from pydantic import BaseModel
from visyn_core import manager


class ReactionCimeSettings(BaseModel):
    dburl: str = "postgresql://admin:admin@localhost:5432/db"
    # statement_timeout: Any = None
    bundles_dir: str | None = None
    # migration: Dict = {"autoUpgrade": True}


# TODO: We can now actually use the type-safe settings...
def get_settings() -> ReactionCimeSettings:
    return manager.settings.reaction_cime  # type: ignore
