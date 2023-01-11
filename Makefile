.DEFAULT_GOAL := help
pkg_src = reaction_cime

black = black --line-length 140 $(pkg_src) setup.py
pyright = pyright $(pkg_src) setup.py
ruff = ruff $(pkg_src) setup.py --line-length 140 --select E,W,F,N,I,C,B,UP,PT,SIM,RUF --ignore E501,C901,B008

.PHONY: start  ## Start the development server
start:
	python $(pkg_src)

.PHONY: all  ## Perform the most common development-time rules
all: format lint test

.PHONY: ci  ## Run all CI validation steps without making any changes to code
ci: check-format lint test

.PHONY: format  ## Auto-format the source code
format:
	$(ruff) --fix
	$(black)

.PHONY: check-format  ## Check the source code format without changes
check-format:
	$(black) --check

.PHONY: lint  ## Run flake8 and pyright
lint:
	$(ruff) --format=github
	$(pyright)

.PHONY: test  ## Run tests
test:
	pytest $(pkg_src)

.PHONEY: documentation ## Generate docs
documentation:
	echo "TODO"

.PHONY: install  ## Install the requirements
install:
	pip install -e .

.PHONY: develop  ## Set up the development environment
develop:
	pip install -e .[develop]

.PHONY: env_encrypt ## Encrypts the current ./<app>/.env
env_encrypt:
	openssl aes-256-cbc -pbkdf2 -in ./$(pkg_src)/.env -out ./$(pkg_src)/.env.enc

.PHONY: env_decrypt ## Decrypts the ./<app>/.env.enc
env_decrypt:
	@if [ -z "${ENV_PASSWORD}" ]; then \
		echo "No ENV_PASSWORD set, prompting for password..."; \
		openssl aes-256-cbc -pbkdf2 -d -in ./$(pkg_src)/.env.enc -out ./$(pkg_src)/.env; \
	else \
		echo "ENV_PASSWORD set, using it..."; \
		openssl aes-256-cbc -pbkdf2 -d -in ./$(pkg_src)/.env.enc -out ./$(pkg_src)/.env -pass env:ENV_PASSWORD; \
	fi

.PHONY: build  ## Build a wheel
build:
	python setup.py sdist bdist_wheel --dist-dir dist_python

.PHONY: publish  ## Publish the ./dist/* using twine
publish:
	pip install twine==3.8.0
	twine upload --repository-url https://upload.pypi.org/legacy/ dist_python/*

.PHONY: help  ## Display this message
help:
	@grep -E \
		'^.PHONY: .*?## .*$$' $(MAKEFILE_LIST) | \
		sort | \
		awk 'BEGIN {FS = ".PHONY: |## "}; {printf "\033[36m%-20s\033[0m %s\n", $$2, $$3}'
