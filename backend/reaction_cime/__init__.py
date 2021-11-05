from .reaction_cime_api import reaction_cime_api
from projection_space_explorer import pse_api


def create_app():
    from flask import Flask
    from flask_cors import CORS
    import logging

    app = Flask(__name__)
    CORS(app)
    app.logger.setLevel(logging.INFO)

    app.register_blueprint(reaction_cime_api)
    return app
