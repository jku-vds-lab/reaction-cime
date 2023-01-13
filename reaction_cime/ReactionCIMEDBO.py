import logging
import time

import pandas as pd
from sqlalchemy import MetaData, Table

_log = logging.getLogger(__name__)


class ReactionCIMEDBO:
    def __init__(self, db):
        self.db = db
        self.metadata = MetaData(db.engine)  # , reflect=True

    def get_table_names(self):
        # both versions work; however, it seems like this delivers the most up-to-date version
        tbl_lst = self.db.engine.table_names()
        if "temp" in tbl_lst:  # we don't wanna show the temp table
            tbl_lst.remove("temp")
        return tbl_lst
        # return self.metadata.tables.keys()

    def save_dataframe(self, df, table_name):
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

        start_time = time.time()
        res = df.to_sql(
            table_name, con=self.db.engine, if_exists="replace", index=True, index_label="id", chunksize=5000
        )  # , schema=table.sql_schema()) # method="multi" 'STRING PRIMARY KEY'
        # https://www.sqlitetutorial.net/sqlite-index/
        delta_time = time.time() - start_time
        _log.info("--- took %i min %f s to save file %s" % (delta_time / 60, delta_time % 60, table_name))

        start_time = time.time()
        self.db.session.execute("DROP INDEX IF EXISTS %s_index;" % table_name)
        delta_time = time.time() - start_time
        _log.info("--- took %i min %f s to drop index of file %s" % (delta_time / 60, delta_time % 60, table_name))

        start_time = time.time()
        self.db.session.execute("create unique index %s_index on %s(id)" % (table_name, table_name))
        delta_time = time.time() - start_time
        _log.info("--- took %i min %f s to create index of file %s" % (delta_time / 60, delta_time % 60, table_name))

        return res

    def update_row(self, table, row_id, value_dict):
        res = self.db.session.query(table).filter(table.c.id == row_id).update(value_dict, synchronize_session=False)
        return res

    # This version of bulk update takes ~15s for 5000 rows
    def update_row_list(self, table_name, id_list, update_dict_list):  # TODO: generalize
        # self.metadata.reflect(only=[table_name])
        table = Table(table_name, self.metadata, autoload=True, autoload_with=self.db.engine)
        for i in range(len(id_list)):
            id = id_list[i]
            update_dict = {}
            for key in update_dict_list:
                update_dict[key] = update_dict_list[key][i]

            self.update_row(table, id, update_dict)

        self.db.session.commit()

    # This version of bulk update is much faster (5000rows: <1s)
    def update_row_bulk(self, table_name, id_list, update_dict_list):
        # https://stackoverflow.com/questions/41870323/sqlalchemy-bulk-update-strategies/41884725#41884725
        # https://stackoverflow.com/questions/19270259/update-with-join-in-sqlite
        df = pd.DataFrame(update_dict_list)
        df.index = id_list
        self.save_dataframe(df, "temp")
        for key in update_dict_list:
            # self.db.session.execute("update %s join temp using(id) set %s.%s = temp.%s"%(table_name,table_name, key, key))
            self.db.session.execute(
                'update "%s" set "%s" = (SELECT "%s" from temp where id = "%s".id)' % (table_name, key, key, table_name)
            )

        self.db.session.commit()

    def get_dataframe_from_table(self, table_name, columns=None):
        if columns is not None and "id" not in columns:
            columns.append("id")
        df = pd.read_sql(table_name, self.db.engine, index_col="id", columns=columns or [])
        if columns is not None:
            columns.remove("id")
            df = df[columns]
        return df

    def get_dataframe_from_table_complete_filter(self, table_name, filter):
        """
        Routes an SQL query including the filter on the given table_name and returns a corresponding pandas dataframe.
        The database in use is fixed as self.db.engine
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
        sql_stmt = 'SELECT * FROM "' + table_name + '" ' + filter
        return pd.read_sql(sql_stmt, self.db.engine, index_col="id")

    def get_dataframe_from_table_filter(self, table_name, filter, columns=None, max_datapoints=-1):
        """
        Routes an SQL query including the filter on the given table_name and returns a corresponding pandas dataframe.
        The database in use is fixed as self.db.engine

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

        sql_stmt = "SELECT " + select_cols + ' FROM "' + table_name + '"'
        if filter is not None and filter != "":
            sql_stmt += " WHERE " + filter
            if max_datapoints > 0:
                datapoint_count = self.get_filter_mask(table_name, filter)["mask"].sum()
                if datapoint_count > max_datapoints:
                    subsample_flag = True
                    # https://www.sqlitetutorial.net/sqlite-functions/sqlite-random/
                    sql_stmt += f" AND abs(RANDOM()%100) < {max_datapoints/datapoint_count*100}"  # add random selection filter based on the ratio between maximum allowed datapoints and NO datapoints that would be selected without sampling

        return pd.read_sql(sql_stmt, self.db.engine, index_col="id"), subsample_flag

    def get_filter_mask(self, table_name, filter):
        sql_stmt = "SELECT id, CASE WHEN " + filter + ' THEN true ELSE false END as mask FROM "' + table_name + '"'
        mask = pd.read_sql(sql_stmt, self.db.engine, index_col="id")
        mask["mask"] = mask["mask"].astype("bool")
        return mask

    def get_value_range_from_table(self, table_name, col_name):
        sql_stmt = f'SELECT min("{col_name}") as min, max("{col_name}") as max FROM "{table_name}"'
        return pd.read_sql(sql_stmt, self.db.engine)

    def get_category_count(self, table_name, col_name):
        sql_stmt = f'SELECT "{col_name}", COUNT(*) as count FROM "{table_name}" GROUP BY "{col_name}"'
        return pd.read_sql(sql_stmt, self.db.engine)

    def get_no_points_from_table(self, table_name):
        sql_stmt = f'SELECT COUNT(*) as count FROM "{table_name}"'
        return pd.read_sql(sql_stmt, self.db.engine)

    def drop_table(self, table_name):
        _log.info("--------drop_table")
        # base = declarative_base()
        self.metadata.reflect(only=[table_name])
        table = self.metadata.tables.get(table_name)
        # _log.info(base)
        _log.info(table)
        # _log.info(base.metadata.tables)
        if table is not None:
            _log.info(f"Deleting {table_name} table")
            # base.metadata.drop_all(self.db.engine, [table], checkfirst=True)
            self.metadata.drop_all(self.db.engine, [table], checkfirst=True)
            return True
            # return table.drop(self.db.engine, checkfirst=True)
        return False
