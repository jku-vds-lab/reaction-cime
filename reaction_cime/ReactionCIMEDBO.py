import logging
import time
from datetime import datetime
from pathlib import Path
from uuid import UUID

import numpy as np
import pandas as pd
from sqlalchemy import MetaData, Table, bindparam, desc
from sqlalchemy.engine import Engine
from visyn_core.security import can_read, current_username

from .db import create_session, get_engine
from .models import Project

_log = logging.getLogger(__name__)


def is_valid_uuid(uuid_to_test):
    try:
        to_test = str(uuid_to_test)
        return str(UUID(to_test)) == to_test
    except ValueError:
        return False


class ReactionCIMEDBO:
    def get_table_names(self):
        with create_session() as session:
            return [
                {
                    "name": p.name,
                    "id": p.id,
                    "creator": p.creator,
                    "permissions": p.permissions,
                    "buddies": p.buddies,
                    "group": p.group,
                    "file_status": p.file_status,
                }
                for p in session.query(
                    Project.id, Project.name, Project.creator, Project.permissions, Project.buddies, Project.group, Project.file_status  # type: ignore
                )
                .order_by(desc(Project.creation_date))
                .all()
                if can_read(p)
            ]

    def get_table_name(self, id, with_schema_prefix: bool = True) -> str:
        if not is_valid_uuid(id):
            raise Exception(f"ID {id} is not a valid UUID.")
        return ("cime4r." if with_schema_prefix else "") + f"data_{id}".replace("-", "_")

    def get_project(self, id) -> Project:
        with create_session() as session:
            project = session.get(Project, id)
            if not project or not can_read(project):
                raise Exception(f"Project with id {id} not found")
            return project

    def update_project(self, id, value_dict) -> bool:
        # add a if can_write and so on
        with create_session() as session:
            res = session.query(Project).filter(Project.id == id).update(value_dict)
            session.commit()
            return res == 1

    def save_project(self, save_name: str):
        p: Project

        with create_session() as session:
            p = Project()
            p.name = save_name  # type: ignore
            p.description = ""  # type: ignore
            # Security
            p.buddies = []  # type: ignore
            p.permissions = 744  # TODO: the user itself can write, everyone else can read
            p.creator = current_username()
            p.creation_date = datetime.utcnow()
            p.modifier = None
            p.modification_date = None
            p.file_status = "not_started"

            session.add(p)
            session.commit()
            session.refresh(p)

        return p

    def save_dataframe(self, path: Path, save_name: str, chunksize: int, project_id, cancel_event):
        msg = ""
        with create_session() as session:
            _log.info("retreiving project")
            p = session.get(Project, project_id)

            p.file_status = "Processing_0"  # type: ignore

            session.commit()

            table_name_with_schema = self.get_table_name(p.id, with_schema_prefix=True)
            table_name = self.get_table_name(p.id, with_schema_prefix=False)

            engine: Engine = session.get_bind()
            with engine.begin() as conn:
                start_time = time.time()

                with pd.read_csv(path, chunksize=chunksize) as reader:
                    for chunk_index, chunk in enumerate(reader):
                        if cancel_event.is_set():
                            raise Exception("cancelled dataframe processing")

                        # Sanitize the column names (i.e. disallow % as it breaks column names in postgres)
                        chunk.columns = [c.strip().replace("%", "_") for c in chunk.columns]
                        if "x" not in chunk.columns or "y" not in chunk.columns:
                            chunk["x"] = np.random.uniform(-50, 50, len(chunk))
                            chunk["y"] = np.random.uniform(-50, 50, len(chunk))
                            msg = "Could not find x and y coordinates in the dataset. Randomly initialized x and y coordinates. Project the dataset to get a better visualization."

                        chunk.to_sql(
                            table_name,
                            schema="cime4r",
                            con=conn,
                            if_exists="replace" if chunk_index == 0 else "append",
                            index=True,
                            index_label="id",
                            chunksize=chunksize,
                        )

                        p.file_status = f"Processing_{chunk_index}"
                        session.commit()

                        _log.info("--- saved chunk %i of file %s" % (chunk_index, save_name))

                p.file_status = "done"  # type: ignore

                delta_time = time.time() - start_time
                _log.info("--- took %i min %f s to save file %s" % (delta_time / 60, delta_time % 60, save_name))

                start_time = time.time()
                conn.execute(f"DROP INDEX IF EXISTS {table_name}_index;")
                delta_time = time.time() - start_time
                _log.info("--- took %i min %f s to drop index of file %s" % (delta_time / 60, delta_time % 60, save_name))

                start_time = time.time()
                conn.execute(f"create unique index {table_name}_index on {table_name_with_schema}(id)")
                delta_time = time.time() - start_time
                _log.info("--- took %i min %f s to create index of file %s" % (delta_time / 60, delta_time % 60, save_name))

            session.commit()

            return str(p.id), msg

    # This version of bulk update is much faster (5000rows: <1s)
    def update_row_bulk(self, id, id_list, update_dict_list: dict):
        df = pd.DataFrame(update_dict_list)
        df.index = id_list  # type: ignore
        # Transform payload to [{_id: 0, x: ..., y: ...}, {_id: 1, x: ..., y: ...}, ...]
        mappings = [{"_id": id, **values} for id, values in df.to_dict("index").items()]
        mappings.sort(key=lambda x: x["_id"])  # type: ignore
        with create_session() as session:
            # Bulk execute the update
            table = Table(
                self.get_table_name(id, with_schema_prefix=False),
                MetaData(session.get_bind(), schema="cime4r"),
                autoload=True,
                autoload_with=session.get_bind(),
            )
            stmt = table.update().where(table.c.id == bindparam("_id")).values({key: bindparam(key) for key in update_dict_list})
            session.execute(stmt, mappings)
            session.commit()

    def get_dataframe_from_table(self, id, columns: list[str] | None = None, chunksize: int | None = None):
        if columns is not None and "id" not in columns:
            columns.append("id")
        with get_engine().connect() as conn:  # NOQA: SIM117
            # We need to set stream_results to true as otherwise read_sql_table will fetch the entire dataframe...
            with conn.execution_options(stream_results=True) as bind:
                chunks = pd.read_sql_table(
                    self.get_table_name(id, with_schema_prefix=False),
                    bind,
                    schema="cime4r",
                    index_col="id",
                    columns=columns or [],
                    chunksize=chunksize,
                )

                # Either return a single dataframe of a generator of dataframes
                for df in [chunks] if isinstance(chunks, pd.DataFrame) else chunks:
                    if columns is not None:
                        df = df.loc[:, [col for col in columns if col != "id"]]  # type: ignore
                    yield df

    def get_dataframe_from_table_complete_filter(self, id, filter):
        """
        Routes an SQL query including the filter on the given table_name and returns a corresponding pandas dataframe.
        The query will look like "SELECT * FROM table_name filter".
        This enables queries such as "SELECT * FROM table_name ORDER BY x LIMIT y".

        Parameters
        ----------
        table_name:
            The table to select from as part of the SQLite query
        filter:
            The filter to apply as part of the SQLite query

        Returns
        ----------
        pandas.Dataframe
            A pandas dataframe containing the query result from the database
        """
        sql_stmt = "SELECT * FROM " + self.get_table_name(id) + " " + filter
        with create_session() as session:
            return pd.read_sql(sql_stmt, session.get_bind(), index_col="id")

    def get_dataframe_from_table_filter(self, id, filter, columns=None, ids: list[int] | None = None, max_datapoints=-1):
        """
        Routes an SQL query including the filter on the given table_name and returns a corresponding pandas dataframe.

        If columns is not None, the query will look like "SELECT columns FROM table_name WHERE filter".
        Otherwise it wil be "SELECT * FROM table_name WHERE filter".

        Parameters
        ----------
        table_name:
            The table to select from as part of the SQLite query
        filter:
            The filter to apply as part of the SQLite query
        columns:
            The columns to perform a select on (SELECT columns FROM)

        max_datapoints:
            if > 0: take a random subsample of the dataset, if there are more than max_datapoints included after the filter

        Returns
        ----------
        pandas.Dataframe
            A pandas dataframe containing the query result from the database
        """
        subsample_flag = False
        select_cols = "*"
        if columns is not None:
            select_cols = "id"
            for col in columns:
                select_cols += f', "{col}"'
        # TODO: This is very bad, as it enables SQL injections! https://xkcd.com/327/
        sql_stmt = "SELECT " + select_cols + " FROM " + self.get_table_name(id)

        all_filters = []

        if filter is not None and filter != "":
            all_filters.append(filter)
            if not ids and max_datapoints > 0:
                datapoint_count = self.get_filter_mask(id, filter)["mask"].sum()
                if datapoint_count > max_datapoints:
                    subsample_flag = True
                    # add random selection filter based on the ratio between maximum allowed datapoints and NO datapoints that would be selected without sampling
                    all_filters.append(f"random()*100 < {max_datapoints/datapoint_count*100} LIMIT {max_datapoints}")

        if ids:
            # If ids are passed, we ignore the filter and only select the ids
            all_filters = [f"id IN ({','.join([str(id) for id in ids])})"]

        if all_filters:
            sql_stmt += " WHERE " + " AND ".join(all_filters)

        with create_session() as session:
            return pd.read_sql(sql_stmt, session.get_bind(), index_col="id"), subsample_flag

    def get_filter_mask(self, id, filter):
        # TODO: This is very bad, as it enables SQL injections! https://xkcd.com/327/
        sql_stmt = "SELECT id, CASE WHEN " + filter + ' THEN true ELSE false END as "mask" FROM ' + self.get_table_name(id) + ""
        with create_session() as session:
            mask = pd.read_sql(sql_stmt, session.get_bind(), index_col="id")
            mask["mask"] = mask["mask"].astype("bool")
            return mask

    def get_value_range_from_table(self, id, col_name):
        # TODO: This is very bad, as it enables SQL injections! https://xkcd.com/327/
        sql_stmt = f'SELECT min("{col_name}") as min, max("{col_name}") as max FROM {self.get_table_name(id)}'
        with create_session() as session:
            return pd.read_sql(sql_stmt, session.get_bind())

    def get_category_count(self, id, col_name):
        # TODO: This is very bad, as it enables SQL injections! https://xkcd.com/327/
        sql_stmt = f'SELECT "{col_name}", COUNT(*) as count FROM {self.get_table_name(id)} GROUP BY "{col_name}"'
        with create_session() as session:
            return pd.read_sql(sql_stmt, session.get_bind())

    def get_no_points_from_table(self, id) -> int:
        sql_stmt = f"SELECT COUNT(*) as count FROM {self.get_table_name(id)}"
        with create_session() as session:
            return int(pd.read_sql(sql_stmt, session.get_bind())["count"][0])

    def drop_table(self, id):
        with create_session() as session:
            p = session.get(Project, id)
            if p:
                _log.info(f"Deleting {id} table")
                session.delete(p)
                session.commit()
                try:
                    session.execute(f"DROP TABLE {self.get_table_name(id)}")
                    session.commit()
                except Exception as e:
                    _log.error(f"Could not drop table {self.get_table_name(id)}: {e}")
                return True
        return False
