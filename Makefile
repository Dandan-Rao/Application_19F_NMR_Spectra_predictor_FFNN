install:
	pip install --upgrade pip &&\
		pip install -r requirements.txt

test:
	## python -m pytest -vv test_predictor.py

# Runs the Black code formatter on all Python files
format:
	black *.py

run:
	python app.py

# Check errors, code quality, and style issues. -R, C focuses on more important errors and warnnings, skipping less critical suggestions
lint:
	pylint --disable=R,C *.py

all: install lint format