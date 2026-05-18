# Makefile
# Tomi Jacobs - Teiko Technical Assessment

# Install all dependencies
setup:
	pip install -r requirements.txt

# Run full pipeline: load data then run all analysis
pipeline:
	python load_data.py
	python analyze.py

# Start the interactive dashboard
dashboard:
	streamlit run dashboard.py
