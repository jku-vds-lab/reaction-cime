import os
from typing import Type

from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from pydantic import BaseModel
from tdp_core.plugin.model import AVisynPlugin, RegHelper

from .settings import ReactionCimeSettings


class VisynPlugin(AVisynPlugin):
    def register(self, registry: RegHelper):
        # registry.append("namespace", "cime_api", "reaction_cime.api", {"namespace": "/api/cime"})

        # registry.append("tdp-sql-database-definition", "reaction_cime", "reaction_cime.db_connector", {"configKey": "reaction_cime"})

        # registry.append(
        #     "tdp-sql-database-migration",
        #     "reaction_cime",
        #     "",
        #     {
        #         "scriptLocation": path.join(path.abspath(path.dirname(__file__)), "migration"),
        #         "configKey": "reaction_cime.migration",
        #         "dbKey": "reaction_cime",
        #     },
        # )

        # Add after server started to cleanup/retry pending processing
        # registry.append("after_server_started", "reaction_cime_retry_datasets", "reaction_cime.after_server_started", {})
        pass

    def init_app(self, app: FastAPI):

        flask_app = Flask(
            __name__,
            template_folder="../../build",
            static_folder="../../build/jku-vds-lab/reaction-cime/static",
            static_url_path="/jku-vds-lab/reaction-cime/static",
        )  # , static_folder="../../build/static", template_folder="../../build")
        # CORS(app)

        basedir = os.path.abspath(os.path.dirname(__file__))
        temp_dir = os.path.join(basedir, "../../temp-files/")
        flask_app.config["tempdir"] = temp_dir
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        flask_app.config["SQLALCHEMY_DATABASE_URI"] = (
            "sqlite:///" + temp_dir + "app.sqlite"
        )  # 'sqlite://' + app.config['REACTION_CIME_FILES_DIRECTORY'] + "//app.sqlite"

        db = SQLAlchemy(flask_app)
        # app.config['SQLALCHEMY_ECHO'] = True

        # setupProjectionTable(db)

        from .ReactionCIMEDBO import ReactionCIMEDBO

        flask_app.config["REACTION_CIME_DBO"] = ReactionCIMEDBO(db)

        # app.config['REACTION_CIME_DBO'].save_dataframe(pd.DataFrame([{"col1":"new", "col2": "test2"}]), "MyTest1")
        # print(app.config['REACTION_CIME_DBO'].get_dataframe_from_table("MyTest1"))
        # app.config['REACTION_CIME_DBO'].update_row("MyTest1", 0, {"col1": "asdf"})
        # print(app.config['REACTION_CIME_DBO'].get_dataframe_from_table("MyTest1"))

        from .reaction_cime_api import reaction_cime_api

        flask_app.register_blueprint(reaction_cime_api)

        app.mount("/api/reaction_cime", WSGIMiddleware(flask_app))

    @property
    def setting_class(self) -> Type[BaseModel]:
        return ReactionCimeSettings
