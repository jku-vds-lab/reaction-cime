import logging

from sqlalchemy import or_

from .db import create_session
from .models import Project

_log = logging.getLogger(__name__)


def create():
    def _():
        with create_session() as session:
            unfinished_projects: list[Project] = (
                session.query(Project)
                .filter(or_(Project.file_status.in_(("processing", "not_started")), Project.file_status.ilike("Processing 18")))  # type: ignore
                .all()
            )

            if unfinished_projects:
                _log.info(f"Found {len(unfinished_projects)} unfinished project(s) since server start, setting file_status to error")

                for project in unfinished_projects:
                    project.file_status = "Error occurred when processing"

                    # TODO: Think about retrying the processing if the file still exists
                    # if project.file_id:
                    #     dataset_processing_executor.submit(process_dataset, project.file_id)

            session.commit()

    return _
