import os

from pydantic import BaseModel
from tdp_core import manager


class ReactionCimeSettings(BaseModel):
    # dburl: str = "sqlite:///:memory:"
    # statement_timeout: Any = None
    tmp_dir: str = os.path.join(os.path.abspath(os.path.dirname(__file__)), "./_data/")
    # migration: Dict = {"autoUpgrade": True}


# TODO: We can now actually use the type-safe settings...
def get_settings() -> ReactionCimeSettings:
    return manager.settings.reaction_cime
