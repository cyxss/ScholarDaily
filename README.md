![logo](ScholarDaily.png)

---
# Description
An automatic tool to collect and summarize the new publishment everyday.
* This is a project for researchers to automaticly scratch the newly-released studies and news happened in the fields they are interested in. Based on the scholar APIs and LLM-based agents, we provide a tool to generated a structural report like a daily newspaper.
* As a researcher in AI4Bio, I provide an example for the topics on Machine learning and Omics data in this version.
* Welcome to provide suggestions or join this project to make ScholarDaily better.
---
# Usage
```python
from scholar_daily import ScholarDaily

sd = ScholarDaily()
# update the config.json before run this step.
sd.load_configs()

# You can set some parameters here. 
sd.set_keywords(["Cardiac Aging", "multiomics","Hypergraph", "Pre-trian"])
sd.set_sources(['Arxiv', 'PubMed']) # In current version of ScholarDaily, you can select sources from  'Arxiv', 'PubMed' and 'bioRxiv'

# run the producing step
sd.produce()
# the output html file can be found in ./ScholarDaily_{date}.html
```
---
# ToDo List
- [x] Add tools for search the Google Scholar.
- [x] Add tools to achieve news.
- [x] Add analysis of each paper.
- [ ] Enable multiple groups of keywords.
- [ ] Improve the style of the html file.

