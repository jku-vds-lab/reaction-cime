# CIME4R

This is the repository for the public cime4r library (as discussed in the paper). It builds upon the PSE library found under the following repository https://github.com/jku-vds-lab/projection-space-explorer.

## Development

The repository is split into frontend (`src`, `package.json`, ...) and backend (`reaction_cime`, `Makefile`, `requirements.txt`, ...). Make sure you have the latest yarn version installed (`corepack enable` when using Node 16).

## Installation

```bash
git clone https://github.com/jku-vds-lab/reaction-cime.git
cd reaction-cime
```

### Frontend

First install the required dependencies

```bash
yarn install
```

and launch the webpack-dev-server via

```bash
yarn start
```

### Backend

First, create a new virtual environment for the dependencies

```bash
python -m venv .venv
```

and activate it

```bash
# Ubuntu
source .venv/bin/activate

# Windows (cmd)
.\.venv\Scripts\activate
```

Then install all dependencies (including dev dependencies)

```bash
make develop
```

and finally start the server

```bash
python reaction_cime
```

As an alternative, you can also user Docker to start the backend. The dockerized backend is managed by `docker compose`:

```
docker compose up
```



PSE package: "git+ssh://git@github.com:jku-vds-lab/projection-space-explorer#develop", "portal:../projection-space-explorer"