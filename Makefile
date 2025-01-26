install:
	pip install --upgrade pip &&\
		pip install -r requirements.txt

test:
	## python -m pytest -vv test_predictor.py

format:
	black *.py

run:
	python main.py

lint:
	pylint --disable=R,C *.py

all: install lint