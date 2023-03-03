"""Add project schema

Revision ID: 81b46bc5282d
Revises: None
Create Date: 2023-02-20 08:54:52.64216

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = "81b46bc5282d"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()

    connection.execute(
        """
        CREATE SCHEMA cime4r;

        CREATE TABLE cime4r.project
        (
            "id" UUID PRIMARY KEY,
            "name" TEXT NOT NULL,
            "description" TEXT NOT NULL,
            "file_exceptions" bytea,
            "file_constraints" bytea,
            /* Security fields */
            "creator" TEXT NOT NULL,
            "creation_date" timestamp with time zone NOT NULL,
            "group" TEXT,
            "buddies" TEXT[] NOT NULL DEFAULT '{}',
            "permissions" integer NOT NULL,
            "modifier" TEXT,
            "modification_date" timestamp with time zone
        );
        """
    )


def downgrade():
    connection = op.get_bind()

    connection.execute(
        """
        DROP SCHEMA cime4r CASCADE;
        """
    )
