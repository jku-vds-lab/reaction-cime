from typing import Any, Dict

from pydantic import BaseModel
from tdp_core import manager


class ReactionCimeSettings(BaseModel):
    dburl: str = "sqlite:///:memory:"
    statement_timeout: Any = None
    uploaded_files_path: str = "/tmp/cime_bayer_files"
    migration: Dict = {"autoUpgrade": True}


# TODO: We can now actually use the type-safe settings...
def get_settings() -> ReactionCimeSettings:
    return manager.settings.reaction_cime
