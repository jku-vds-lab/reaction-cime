"""Add fully processed flag

Revision ID: 027d77624694
Revises: 81b46bc5282d
Create Date: 2023-03-31 10:26:52.129803

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = "027d77624694"
down_revision = "81b46bc5282d"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()

    connection.execute(
        """
        ALTER TABLE cime4r.project
        ADD "fully_processed" BOOLEAN NOT NULL DEFAULT false;
        """
    )


def downgrade():
    connection = op.get_bind()

    connection.execute(
        """
        ALTER TABLE cime4r.project
        DROP COLUMN "fully_processed";
        """
    )
