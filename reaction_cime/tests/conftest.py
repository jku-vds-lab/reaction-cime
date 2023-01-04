import os
from typing import Any, Generator

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from tdp_core.security.manager import SecurityManager
from tdp_core.security.model import User
from tdp_core.server.visyn_server import create_visyn_server


@pytest.fixture(scope="session")
def app() -> Generator[FastAPI, Any, None]:
    server = create_visyn_server(
        workspace_config={
            "_env_file": os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.env"),
            "tdp_core": {"enabled_plugins": ["reaction_cime"]},
            "reaction_cime": {},
        }
    )

    yield server


@pytest.fixture()
def client(monkeypatch, app: FastAPI):
    def mock_current_user_in_manager(self):
        return User(id="admin")

    monkeypatch.setattr(SecurityManager, "current_user", property(mock_current_user_in_manager))

    with TestClient(app) as client:
        yield client
