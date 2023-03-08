import logging
import time
from datetime import datetime

import pandas as pd
from sqlalchemy import MetaData, Table, bindparam
from sqlalchemy.engine import Engine
from visyn_core.security import can_read, current_username

from .db import create_session
from .models import Project

_log = logging.getLogger(__name__)


class ReactionCIMEDBO:
    def get_table_names(self):
        with create_session() as session:
            return [{"name": p.name, "id": p.id, "creator": p.creator, "permissions": p.permissions, "buddies": p.buddies, "group": p.group} for p in session.query(Project.id, Project.name, Project.creator, Project.permissions, Project.buddies, Project.group).all() if can_read(p)]  # type: ignore

    def get_table_name(self, id, with_schema_prefix: bool = True) -> str:
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

    def save_dataframe(self, df: pd.DataFrame, filename: str) -> str:
        # primary key is not set to index of dataframe by default. this code could be a workaround, but does not work yet
        # https://stackoverflow.com/questions/30867390/python-pandas-to-sql-how-to-create-a-table-with-a-primary-key/69350803#69350803
        # https://github.com/pandas-dev/pandas/blob/945c9ed766a61c7d2c0a7cbb251b6edebf9cb7d5/pandas/io/sql.py#L933
        # pandas_sql = pd.io.sql.pandasSQL_builder(self.db.engine)
        # table = SQLiteTable(
        #     table_name,
        #     pandas_sql,
        #     frame=df,
        #     index=True,
        #     index_label="id",
        #     keys="id",
        # )
        # # _log.info(table.sql_schema())
        # table.create()
        # res = table.insert(chunksize=5000, method=None) # method="multi"
        with create_session() as session:
            p = Project()
            p.name = filename  # type: ignore
            p.description = ""  # type: ignore
            # Security
            p.buddies = []  # type: ignore
            p.permissions = 744  # TODO: the user itself can write, everyone else can read
            p.creator = current_username()
            p.creation_date = datetime.utcnow()
            p.modifier = None
            p.modification_date = None
            session.add(p)
            session.flush()
            session.refresh(p)

            table_name_with_schema = self.get_table_name(p.id, with_schema_prefix=True)
            table_name = self.get_table_name(p.id, with_schema_prefix=False)

            engine: Engine = session.get_bind()
            with engine.begin() as conn:
                start_time = time.time()
                df.to_sql(
                    table_name, schema="cime4r", con=conn, if_exists="replace", index=True, index_label="id", chunksize=5000
                )  # , schema=table.sql_schema()) # method="multi" 'STRING PRIMARY KEY'
                # https://www.sqlitetutorial.net/sqlite-index/
                delta_time = time.time() - start_time
                _log.info("--- took %i min %f s to save file %s" % (delta_time / 60, delta_time % 60, filename))

                start_time = time.time()
                conn.execute(f"DROP INDEX IF EXISTS {table_name}_index;")
                delta_time = time.time() - start_time
                _log.info("--- took %i min %f s to drop index of file %s" % (delta_time / 60, delta_time % 60, filename))

                start_time = time.time()
                conn.execute(f"create unique index {table_name}_index on {table_name_with_schema}(id)")
                delta_time = time.time() - start_time
                _log.info("--- took %i min %f s to create index of file %s" % (delta_time / 60, delta_time % 60, filename))

            session.commit()

            return str(p.id)

    # This version of bulk update is much faster (5000rows: <1s)
    def update_row_bulk(self, id, id_list, update_dict_list: dict):
        df = pd.DataFrame(update_dict_list)
        df.index = id_list  # type: ignore
        # Transform payload to [{_id: 0, x: ..., y: ...}, {_id: 1, x: ..., y: ...}, ...]
        mappings = [{"_id": id, **values} for id, values in df.to_dict("index").items()]
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

    def get_dataframe_from_table(self, id, columns=None):
        if columns is not None and "id" not in columns:
            columns.append("id")
        with create_session() as session:
            df = pd.read_sql_table(
                self.get_table_name(id, with_schema_prefix=False),
                session.get_bind(),
                schema="cime4r",
                index_col="id",
                columns=columns or [],
            )
            if columns is not None:
                columns.remove("id")
                df = df[columns]
            return df

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

    def get_dataframe_from_table_filter(self, id, filter, columns=None, max_datapoints=-1):
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
        sql_stmt = "SELECT " + select_cols + " FROM " + self.get_table_name(id) + ""
        if filter is not None and filter != "":
            sql_stmt += " WHERE " + filter
            if max_datapoints > 0:
                datapoint_count = self.get_filter_mask(id, filter)["mask"].sum()
                if datapoint_count > max_datapoints:
                    subsample_flag = True
                    # https://www.sqlitetutorial.net/sqlite-functions/sqlite-random/
                    sql_stmt += f" AND abs(RANDOM()%100) < {max_datapoints/datapoint_count*100}"  # add random selection filter based on the ratio between maximum allowed datapoints and NO datapoints that would be selected without sampling
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

    def get_no_points_from_table(self, id):
        # TODO: This is very bad, as it enables SQL injections! https://xkcd.com/327/
        sql_stmt = f"SELECT COUNT(*) as count FROM {self.get_table_name(id)}"
        with create_session() as session:
            return pd.read_sql(sql_stmt, session.get_bind())

    def drop_table(self, id):
        with create_session() as session:
            p = session.get(Project, id)
            if p:
                _log.info(f"Deleting {id} table")
                session.delete(p)
                # TODO: This is very bad, as it enables SQL injections! https://xkcd.com/327/
                session.execute(f"DROP TABLE {self.get_table_name(id)}")
                session.commit()
                return True
        return False
