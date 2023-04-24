"""add_file_status_field

Revision ID: 20b7a0330880
Revises: 027d77624694
Create Date: 2023-04-21 10:37:21.857830

"""
from alembic import op


# revision identifiers, used by Alembic.
revision = '20b7a0330880'
down_revision = '027d77624694'
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()

    connection.execute(
        """
        ALTER TABLE cime4r.project
        ADD "file_status" TEXT;

        UPDATE cime4r.project
        SET "file_status" = 'done';

        ALTER TABLE cime4r.project
        ALTER COLUMN "file_status" SET NOT NULL;

        ALTER TABLE cime4r.project
        DROP COLUMN "fully_processed";
        """
    )
    pass


def downgrade():
    connection = op.get_bind()

    connection.execute(
        """
        ALTER TABLE cime4r.project
        DROP COLUMN "file_status";

        ALTER TABLE cime4r.project
        ADD "fully_processed" BOOLEAN NOT NULL DEFAULT false;
        """
    )
    pass
