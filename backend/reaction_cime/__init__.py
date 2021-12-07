from .reaction_cime_api import reaction_cime_api
from projection_space_explorer import pse_api
from .db import ACimeDBO
from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy
import os
import logging
import pandas as pd
from sqlalchemy import MetaData, Table, update, String, Column, Integer
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import Session


def create_app():
    print("---create_app-----")

    app = Flask(__name__)
    CORS(app)

    basedir = os.path.abspath(os.path.dirname(__file__))
    temp_dir = os.path.join(basedir, 'temp-files/')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + temp_dir + 'app.sqlite' #'sqlite://' + app.config['REACTION_CIME_FILES_DIRECTORY'] + "//app.sqlite"

    app.logger.setLevel(logging.INFO)

    db = SQLAlchemy(app)
    db.create_all()

    app.config['REACTION_CIME_DBO'] = ReactionCIMEDBO(db)

    # app.config['REACTION_CIME_DBO'].save_dataframe(pd.DataFrame([{"col1":"new", "col2": "test2"}]), "MyTest1")
    # print(app.config['REACTION_CIME_DBO'].get_dataframe_from_table("MyTest1"))
    # app.config['REACTION_CIME_DBO'].update_row("MyTest1", 0, {"col1": "asdf"})
    # print(app.config['REACTION_CIME_DBO'].get_dataframe_from_table("MyTest1"))

    print("---")

    app.register_blueprint(reaction_cime_api)
    return app


class ReactionCIMEDBO():

    def __init__(self, db):
        self.db = db
        self.metadata = MetaData(db.engine, reflect=True)

    def get_table_names(self):
        # both versions work; however, it seems like this delivers the most up-to-date version
        return self.db.engine.table_names()
        # return self.metadata.tables.keys()

    def save_dataframe(self, df, table_name):
        # to_sql(self.db.engine, df, table_name, if_exists='replace')
        return df.to_sql(table_name, con=self.db.engine, if_exists="replace", index=True, chunksize=5000)

    def update_row(self, table_name, row_id, value_dict):        
        # self.metadata.reflect(only=[table_name])
        table = Table(table_name, self.metadata, autoload=True, autoload_with=self.db.engine)
        
        res = self.db.session.query(table).filter(table.c.index == row_id).update(value_dict, synchronize_session = False)
        self.db.session.commit()

        return res

    def get_dataframe_from_table(self, table_name, columns=None):
        if columns is not None and "index" not in columns:
            columns.append("index")
        return pd.read_sql(table_name, self.db.engine, index_col="index", columns=columns)

    def get_dataframe_from_table_filter(self, table_name, filter):
        sql_stmt = "SELECT * FROM " + table_name + " WHERE " + filter
        return pd.read_sql(sql_stmt, self.db.engine, index_col="index")
    
    def drop_table(self, table_name):
        base = declarative_base()
        table = self.metadata.tables.get(table_name)
        if table is not None:
            logging.info(f'Deleting {table_name} table')
            base.metadata.drop_all(self.db.engine, [table], checkfirst=True)
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