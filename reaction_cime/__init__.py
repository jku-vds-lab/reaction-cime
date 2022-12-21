import os
from typing import Type

from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware
from fastapi.staticfiles import StaticFiles
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from pydantic import BaseModel
from tdp_core.plugin.model import AVisynPlugin, RegHelper
from tdp_core.server.utils import init_legacy_app

from .settings import ReactionCimeSettings, get_settings


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
        import logging

        _log = logging.getLogger(__name__)

        flask_app = Flask(__name__)
        # CORS(app)

        tmp_dir = get_settings().tmp_dir
        flask_app.config["tempdir"] = tmp_dir
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        flask_app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///" + os.path.join(tmp_dir, "app.sqlite")

        db = SQLAlchemy(flask_app)
        # app.config['SQLALCHEMY_ECHO'] = True

        # setupProjectionTable(db)

        from .ReactionCIMEDBO import ReactionCIMEDBO

        flask_app.config["REACTION_CIME_DBO"] = ReactionCIMEDBO(db)

        from .reaction_cime_api import reaction_cime_api

        flask_app.register_blueprint(reaction_cime_api)

        init_legacy_app(flask_app)
        app.mount("/api/reaction_cime", WSGIMiddleware(flask_app))

        @app.on_event("startup")
        async def startup():
            # Add the / path at the very end to match all other routes before
            bundles_dir = get_settings().bundles_dir
            if bundles_dir:
                # Mount the bundles directory as static file to enable the frontend (required in single Dockerfile mode)
                _log.info(f"Mounting bundles dir: {bundles_dir}")
                app.mount("/", StaticFiles(directory=bundles_dir, html=True), name="reaction_cime_bundles")

    @property
    def setting_class(self) -> Type[BaseModel]:
        return ReactionCimeSettings
