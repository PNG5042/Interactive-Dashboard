.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python initial_analysis.py
	python statistical_analysis.py
	python data_subset_analysis.py

dashboard:
	streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501
