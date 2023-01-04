from abc import abstractmethod

from sqlalchemy.orm import Session, scoped_session
from tdp_core import manager


class ACimeDBO:
    def get_datasets(self):
        return self.get_datasets_by()

    @abstractmethod
    def get_datasets_by(self, **kwargs):
        pass

    @abstractmethod
    def get_dataset_by(self, **kwargs):
        pass

    @abstractmethod
    def save_dataset(self, dataset):
        pass

    @abstractmethod
    def delete_dataset_by(self, **kwargs):
        pass

    @abstractmethod
    def save_file(self, file):
        pass


def get_engine():
    return manager.db.engine("dv_reaction_cime")


def create_web_session() -> Session:
    return manager.db.create_web_session(get_engine())


def create_session() -> Session:
    return manager.db.create_session(get_engine())


def create_scoped_session() -> Session:
    return scoped_session(manager.db._sessionmakers[get_engine()])()
