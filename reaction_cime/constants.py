from pathlib import Path

from .settings import get_settings

storage_path: Path = Path(get_settings().uploaded_files_path) / "storage"
