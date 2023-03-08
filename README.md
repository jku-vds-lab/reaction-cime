# CIME4R

This is the repository for the public cime4r library (as discussed in the paper). It builds upon the PSE library found under the following repository https://github.com/jku-vds-lab/projection-space-explorer.

## Development

The repository is split into frontend (`src`, `package.json`, ...) and backend (`reaction_cime`, `Makefile`, `requirements.txt`, ...). Make sure you have the latest yarn version installed (`corepack enable` when using Node 16).

## Installation

```bash
git clone https://github.com/jku-vds-lab/reaction-cime.git
cd reaction-cime
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

start the database

```
docker compose up
```

and finally start the server

```bash
python reaction_cime
```


Now, if a login screen pops up, you can use admin:admin to login. If you want to disable the login screen and go directly to the application, create a `reaction_cime/.env` with the following contents. After restaring the server, you will be automatically logged in.

```
VISYN_CORE__SECURITY__STORE__NO_SECURITY_STORE__ENABLE=true
VISYN_CORE__SECURITY__STORE__NO_SECURITY_STORE__USER=admin
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


### Deploy with Docker

```bash
yarn install
```

```bash
yarn run webpack:prod
```

```bash
docker build -f Dockerfile -t reaction_cime .
```

```bash
docker run --rm -it --network host reaction_cime
```

Beware that you will need a Postgres to run the image. By default, it will use the connection string in `settings.py`, which you can override via ENV variables. For example, you can set `REACTION_CIME__DBURL=postgresql://...` and use any database of your liking.

## Linking PSE and Reaction-CIME Frontend

If you want to make changes to PSE and view the changes without having to push to the repo and reinstalling dependencies, the recommended way is to use the yarn link/portal and/or our webpack resolveAliases feature.

First, clone `projection-space-explorer` into the current directory (i.e. into the reaction-cime directory). Do not install `projection-space-explorer`, as we do not want any `node_modules` within that folder, as it should use the ones from the reaction-cime directory.

Add a portal to the local `projection-space-explorer` in the reaction-cime package.json (**this is a local change and should not be committed!**):

```json
  "resolutions": {
    ...
    "projection-space-explorer": "portal:./projection-space-explorer"
  },
```

Now, install everything via `yarn install`.

To now include `projection-space-explorer` to your webpack build, add a `.yo-rc-workspace.json` and update the `resolveAliases` (note the single `.`, i.e. using the current folder):

```json
{
  "resolveAliases": {
    "projection-space-explorer": "./projection-space-explorer/src/index.ts"
  }
}
```

With that, you can now edit all files of `projection-space-explorer`, including auto-completion (as the node_modules of the application will be used as main lookup), and get hot-reloading.
