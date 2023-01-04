import os

from tdp_core.server.visyn_server import create_visyn_server

app = create_visyn_server(workspace_config={"_env_file": os.path.join(os.path.dirname(os.path.realpath(__file__)), ".env")})
