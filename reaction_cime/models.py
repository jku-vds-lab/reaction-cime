import logging
import pickle
import uuid
import zlib
from datetime import datetime
from uuid import UUID as UUID_TYPE

from sqlalchemy import ARRAY, TEXT, Boolean, Column, DateTime, Integer, PickleType, String
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

_log = logging.getLogger(__name__)


class CompressedPickler:
    def loads(self, data: bytes, **kwargs) -> object:
        return pickle.loads(zlib.decompress(data))

    def dumps(self, data: bytes, protocol=None, **kwargs) -> bytes:
        return zlib.compress(pickle.dumps(data, protocol=protocol), level=zlib.Z_BEST_COMPRESSION)


class Project(Base):  # type: ignore
    __tablename__ = "project"
    __table_args__ = {"schema": "cime4r", "extend_existing": True}

    id: UUID_TYPE = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)  # type: ignore
    name = Column(String, nullable=False)
    description = Column(String, nullable=False, default="")
    file_exceptions = Column(PickleType(pickler=CompressedPickler()))  # type: ignore
    file_constraints = Column(PickleType(pickler=CompressedPickler()))  # type: ignore
    fully_processed = Column(Boolean, nullable=False, default=False)  # type: ignore
    # Security
    creator: str = Column(TEXT, nullable=False)  # type: ignore
    creation_date: datetime = Column(DateTime, nullable=False)  # type: ignore
    group: str | None = Column(TEXT)  # type: ignore
    permissions: int = Column(Integer, nullable=False)  # type: ignore
    buddies = Column(ARRAY(TEXT), nullable=False, default=[])
    modifier: str | None = Column(TEXT)  # type: ignore
    modification_date: datetime | None = Column(DateTime)  # type: ignore
