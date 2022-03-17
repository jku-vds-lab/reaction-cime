from .reaction_cime_api import reaction_cime_api
from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy
import os
import logging
import pandas as pd
from sqlalchemy import MetaData, Table
# from sqlalchemy.ext.declarative import declarative_base


def create_app():
    print("---create_app-----")

    app = Flask(__name__, template_folder="../../build", static_folder="../../build/jku-vds-lab/reaction-cime/static", static_url_path="/jku-vds-lab/reaction-cime/static")#, static_folder="../../build/static", template_folder="../../build")
    CORS(app)

    basedir = os.path.abspath(os.path.dirname(__file__))
    temp_dir = os.path.join(basedir, '../../temp-files/')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + temp_dir + 'app.sqlite' #'sqlite://' + app.config['REACTION_CIME_FILES_DIRECTORY'] + "//app.sqlite"
    # app.config['SQLALCHEMY_ECHO'] = True

    app.logger.setLevel(logging.INFO) # DEBUG

    db = SQLAlchemy(app)
    db.create_all()

    # setupProjectionTable(db)

    app.config['REACTION_CIME_DBO'] = ReactionCIMEDBO(db)

    # app.config['REACTION_CIME_DBO'].save_dataframe(pd.DataFrame([{"col1":"new", "col2": "test2"}]), "MyTest1")
    # print(app.config['REACTION_CIME_DBO'].get_dataframe_from_table("MyTest1"))
    # app.config['REACTION_CIME_DBO'].update_row("MyTest1", 0, {"col1": "asdf"})
    # print(app.config['REACTION_CIME_DBO'].get_dataframe_from_table("MyTest1"))

    print("---")

    app.register_blueprint(reaction_cime_api)
    return app

# def setupProjectionTable(db):
#     class ProjectionTable(db.Model):
#         __tablename__ = "projection"

#         id = Column(Integer, primary_key=True)
#         projection_name = Column(String, nullable=False)
#         dataset_name = Column(String, nullable=False)
#         iteration_step = Column(Integer, nullable=False)
#         projection_coords = Column(PickleType, nullable=False)


class ReactionCIMEDBO():

    def __init__(self, db):
        self.db = db
        self.metadata = MetaData(db.engine) #, reflect=True

    def get_table_names(self):
        # both versions work; however, it seems like this delivers the most up-to-date version
        tbl_lst = self.db.engine.table_names()
        if "temp" in tbl_lst: # we don't wanna show the temp table
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
        # # print(table.sql_schema())
        # table.create()
        # res = table.insert(chunksize=5000, method=None) # method="multi"
        
        res = df.to_sql(table_name, con=self.db.engine, if_exists="replace", index=True, index_label="id", chunksize=5000)#, schema=table.sql_schema()) # method="multi" 'STRING PRIMARY KEY'
        # https://www.sqlitetutorial.net/sqlite-index/
        self.db.session.execute('DROP INDEX IF EXISTS %s_index;'%table_name)
        self.db.session.execute('create unique index %s_index on %s(id)'%(table_name, table_name))
        return res

    def update_row(self, table, row_id, value_dict):
        res = self.db.session.query(table).filter(table.c.id == row_id).update(value_dict, synchronize_session = False)
        return res

    # This version of bulk update takes ~15s for 5000 rows
    def update_row_list(self, table_name, id_list, update_dict_list): # TODO: generalize
        # self.metadata.reflect(only=[table_name])
        table = Table(table_name, self.metadata, autoload=True, autoload_with=self.db.engine)
        for i in range(len(id_list)):
            id = id_list[i]
            update_dict = {}
            for key in update_dict_list.keys():
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
        for key in update_dict_list.keys():
            # self.db.session.execute("update %s join temp using(id) set %s.%s = temp.%s"%(table_name,table_name, key, key))
            self.db.session.execute('update "%s" set "%s" = (SELECT "%s" from temp where id = "%s".id)'%(table_name, key, key, table_name))

        self.db.session.commit()

    def get_dataframe_from_table(self, table_name, columns=None):
        if columns is not None and "id" not in columns:
            columns.append("id")
        df = pd.read_sql(table_name, self.db.engine, index_col="id", columns=columns)
        if(columns is not None):
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

    def get_dataframe_from_table_filter(self, table_name, filter, columns=None):
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


        Returns
        ----------
        pandas.Dataframe
            A pandas dataframe containing the query result from the database
        """
        select_cols = "*"
        if columns is not None:
            select_cols = "id"
            for col in columns:
                select_cols += f', "{col}"'

        sql_stmt = 'SELECT ' + select_cols + ' FROM "' + table_name + '" WHERE ' + filter
        return pd.read_sql(sql_stmt, self.db.engine, index_col="id")

    def get_filter_mask(self, table_name, filter):
        sql_stmt = 'SELECT id, CASE WHEN ' + filter + ' THEN true ELSE false END as mask FROM "' + table_name + '"'
        mask = pd.read_sql(sql_stmt, self.db.engine, index_col="id")
        mask["mask"] = mask["mask"].astype("bool")
        return mask

    def get_value_range_from_table(self, table_name, col_name):
        sql_stmt = f'SELECT min("{col_name}") as min, max("{col_name}") as max FROM "{table_name}"'
        return pd.read_sql(sql_stmt, self.db.engine)

    def drop_table(self, table_name):
        print("--------drop_table")
        # base = declarative_base()
        self.metadata.reflect(only=[table_name])
        table = self.metadata.tables.get(table_name)
        # print(base)
        print(table)
        # print(base.metadata.tables)
        if table is not None:
            logging.info(f'Deleting {table_name} table')
            # base.metadata.drop_all(self.db.engine, [table], checkfirst=True)
            self.metadata.drop_all(self.db.engine, [table], checkfirst=True)
            return True
            # return table.drop(self.db.engine, checkfirst=True)
        return False

import io
def to_sql(engine, df, table, if_exists='fail', sep='\t', encoding='utf8'):
    # Create Table
    df[:0].to_sql(table, engine, if_exists=if_exists)

    # Prepare data
    output = io.StringIO()
    df.to_csv(output, sep=sep, header=False, encoding=encoding)
    output.seek(0)

    # Insert data
    connection = engine.raw_connection()
    cursor = connection.cursor()
    cursor.copy_from(output, table, sep=sep, null='') # only postgreSQL has this "copy_from" function
    connection.commit()
    cursor.close()