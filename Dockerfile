
# define our environment
FROM python:3.7-buster

# # copy the pre-built front-end --> comment for development because we mount the volume anyway
# COPY build/ /app/build/jku-vds-lab/reaction-cime

# copy everything from our backend to our app folder # need to copy backend because we have to install the python packages
COPY backend/ /app/backend/

# copy everything from temp-files folder that includes the database --> comment for development because we mount the volume anyway
# COPY temp-files/ /app/temp-files/

# define target folder
WORKDIR /app/backend


RUN pip install -e .
ENV FLASK_APP reaction_cime
ENV FLASK_ENV development



CMD [ "python3", "-m", "flask", "run", "--host=0.0.0.0", "--eager-loading"] 

# Running
# docker build -f Dockerfile -t reaction_cime .
# docker run -d -p 5000:5000 --detach reaction_cime
# use this for development (file sharing); first -v contains the "backend" folder; second -v contains the "temp-files" folder where our database is; third -v contains the "build" folder where the front-end gets compiled to
# docker run -d -p 5000:5000 -v "C:/Users/Christina/Workspace/ICG/bayer/reaction-cime/backend:/app/backend" -v "C:/Users/Christina/Workspace/ICG/bayer/reaction-cime/temp-files:/app/temp-files" -v "C:/Users/Christina/Workspace/ICG/bayer/reaction-cime/build:/app/build/jku-vds-lab/reaction-cime"  --detach reaction_cime
