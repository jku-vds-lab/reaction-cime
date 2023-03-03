# This Dockerfile is used as standalone container for simple deployments, it will be built and pushed to https://github.com/orgs/jku-vds-lab/packages?repo_name=reaction-cime automatically by GH Actions in the build.yml
FROM python:3.10-buster

# copy everything from our backend to our app folder # need to copy backend because we have to install the python packages
COPY reaction_cime/ /app/reaction_cime/
COPY Makefile MANIFEST.in README.md setup.py setup.cfg package.json requirements.txt requirements_dev.txt /app/

# define target folder
WORKDIR /app/

# Install some build tools and finally python dependencies (numpy is required to build opentsne)
RUN pip install numpy && make install

# Override the setttings.py to use include the bundled frontend
ENV REACTION_CIME__BUNDLES_DIR /app/bundles
# Disable the login and always use a anonymous user
ENV VISYN_CORE__SECURITY__STORE__NO_SECURITY_STORE__ENABLE true
ENV VISYN_CORE__SECURITY__STORE__NO_SECURITY_STORE__USER anonymous

# copy the pre-built front-end --> comment for development because we mount the volume anyway
COPY bundles/ /app/bundles/

# expose default port
EXPOSE 9000

CMD ["uvicorn", "visyn_core.server.main:app", "--host", "0.0.0.0", "--port", "9000"]

# Running
# docker build -f Dockerfile -t reaction_cime .
# docker run --rm -it -p 9000:9000 reaction_cime
# use this for development (file sharing); first -v contains the "backend" folder; second -v contains the "build" folder where the front-end gets compiled to
# docker run --rm -it -p 9000:9000 -v "$PWD/reaction_cime/:/app/reaction_cime/" -v "$PWD/bundles/:/app/bundles/" reaction_cime