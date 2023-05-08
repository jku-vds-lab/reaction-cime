from os import path

from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware
from fastapi.staticfiles import StaticFiles
from flask import Flask
from pydantic import BaseModel
from visyn_core.plugin.model import AVisynPlugin, RegHelper
from visyn_core.server.utils import init_legacy_app
from visyn_core.settings.client_config import visyn_client_config

from .settings import ReactionCimeSettings, get_settings


class VisynPlugin(AVisynPlugin):
    def register(self, registry: RegHelper):
        # registry.append("namespace", "cime_api", "reaction_cime.api", {"namespace": "/api/cime"})

        registry.append("tdp-sql-database-definition", "reaction_cime", "", {"configKey": "reaction_cime"})

        registry.append(
            "tdp-sql-database-migration",
            "reaction_cime",
            "",
            {
                "scriptLocation": path.join(path.abspath(path.dirname(__file__)), "migration"),
                "configKey": "reaction_cime.migration",
                "dbKey": "reaction_cime",
            },
        )

        # Add after server started to cleanup/retry pending processing
        registry.append("after_server_started", "reaction_cime_retry_datasets", "reaction_cime.after_server_started", {})
        pass

    def init_app(self, app: FastAPI):
        import logging

        _log = logging.getLogger(__name__)

        get_settings().uploaded_files_path.mkdir(parents=True, exist_ok=True)

        flask_app = Flask(__name__)

        from .ReactionCIMEDBO import ReactionCIMEDBO

        dbo = ReactionCIMEDBO()

        # Make DBO accessible in FastAPI
        app.state.reaction_cime_dbo = dbo

        from .reaction_cime_api import cancel_events, project_executor, reaction_cime_api, router

        app.include_router(router, prefix="/api/reaction_cime/v2")

        flask_app.register_blueprint(reaction_cime_api)

        init_legacy_app(flask_app)
        app.mount("/api/reaction_cime", WSGIMiddleware(flask_app))

        @visyn_client_config
        class CIME4RClientConfigModel(BaseModel):
            publicVersion: bool = False  # NOQA N815

        @app.on_event("startup")
        async def startup():
            # Add the / path at the very end to match all other routes before
            bundles_dir = get_settings().bundles_dir
            if bundles_dir:
                # Mount the bundles directory as static file to enable the frontend (required in single Dockerfile mode)
                _log.info(f"Mounting bundles dir: {bundles_dir}")
                app.mount("/", StaticFiles(directory=bundles_dir, html=True), name="reaction_cime_bundles")

        @app.on_event("shutdown")
        async def shutdown():
            for event in cancel_events:
                event.set()
            project_executor.shutdown(wait=True)

    @property
    def setting_class(self) -> type[BaseModel]:
        return ReactionCimeSettings
