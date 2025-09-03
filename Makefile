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
# Adding || true makes the command always exit with zero, so the job continues.
lint:
	pylint --disable=R,C *.py || true

all: install lint format

# Lambda deployment testing targets
test-lambda-quick:
	@echo "🧪 Running quick Lambda deployment test..."
	python quick_lambda_test.py

test-lambda-full:
	@echo "🧪 Running full Lambda deployment test..."
	python test_lambda_deployment.py

test-lambda: test-lambda-quick

# Build and test Lambda deployment
build-and-test-lambda:
	@echo "🔨 Building Docker image..."
	docker build -t nmr-predictor:latest .
	@echo "🧪 Testing Lambda deployment..."
	python quick_lambda_test.py

# Clean up test containers
clean-lambda-test:
	@echo "🧹 Cleaning up Lambda test containers..."
	docker stop lambda-quick-test 2>/dev/null || true
	docker rm lambda-quick-test 2>/dev/null || true
	docker stop lambda-test-container 2>/dev/null || true
	docker rm lambda-test-container 2>/dev/null || true