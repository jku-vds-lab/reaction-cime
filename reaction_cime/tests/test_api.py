from fastapi.testclient import TestClient


def test_loggedinas(client: TestClient):
    assert client.get("/loggedinas").json()["name"] == "admin"
