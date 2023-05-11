from pathlib import Path

from pydantic import BaseModel
from visyn_core import manager


class ReactionCimeSettings(BaseModel):
    dburl: str = "postgresql://admin:admin@localhost:5432/db"
    # statement_timeout: Any = None
    bundles_dir: str | None = None
    uploaded_files_path: Path = Path("/tmp/uploaded_files")
    max_memory_usage_for_projections: int = 50
    """
    Limit of dataset size in MB which we use to compute chunk sizes.
    The real memory used is higher than that, as we allocate new arrays, ...
    """
    # migration: Dict = {"autoUpgrade": True}


# TODO: We can now actually use the type-safe settings...
def get_settings() -> ReactionCimeSettings:
    return manager.settings.reaction_cime  # type: ignore
