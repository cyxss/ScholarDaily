import json
import os
from typing import List
from datetime import date, timedelta, datetime
import arxiv
from tqdm import tqdm
import requests
from bs4 import BeautifulSoup
import time
from metapub import PubMedFetcher
from Bio import Entrez
import markdown
from openai import OpenAI
import pickle

class ScholarDaily:
    def __init__(self, json_path='scholar_daily/config.json'):
        self.json_path = json_path
        self.configs = {}
        self.date = date.today() - timedelta(days = 1)#看昨天一天的
        self.source_results = {}        
        self.All_result = []
        self.fuction_mapping = {'BioRxiv':self.search_bioRxiv,
                                'Arxiv': self.search_Arxiv,
                                'PubMed':self.search_PubMed,
                               "GoogleScholar":self.search_GoogleScholar}
        
    ########## functions on configs ##########
    def load_configs(self):
        """读取JSON文件并将其内容加载为内部变量"""
        if os.path.exists(self.json_path):
            with open(self.json_path, 'r', encoding='utf-8') as f:
                self.configs = json.load(f)
        else:
            raise FileNotFoundError(f"配置文件 {self.json_path} 不存在")

    def save_configs(self):
        """将内部变量保存到JSON文件"""
        with open(self.json_path, 'w', encoding='utf-8') as f:
            json.dump(self.configs, f, ensure_ascii=False, indent=4)

    def set_keywords(self, keywords:List[str]):
        """设置JSON中的keywords"""
        self.configs["Keywords"] = keywords
        self.save_configs()
        
    def set_sources(self, sources:List[str]):
        """设置JSON中的resources"""
        self.configs["Sources"] = sources
        self.save_configs()

    def get_all_configs(self):
        """获取所有参数"""
        return self.configs
    
    def set_date(self, date_str:str):
        try:
            self.date = datetime.strptime(date_str, '%Y%m%d').date()
        except Exception as e:
            print(e)
    ##############################
    
    
    ########## functions about query resource ##########
    def search_GoogleScholar(keywords_list:List[str], target_date):
        print("Coming soon.")

    def search_Arxiv(self):
        keywords_list = self.configs['Keywords']
        target_date = self.date
        
        # Construct the default API client.
        client = arxiv.Client()
        
        # Construct query condition
        query = " OR ".join([f"all:{key}" for key in keywords_list])
        print(f"query condition as:\n {query}")
        
        # Search for the 10 most recent articles matching the keyword "quantum."
        search = arxiv.Search(
          query = query,
          max_results = self.configs['max_num_per_source'],
          sort_by = arxiv.SortCriterion.SubmittedDate
        )

        query_results = client.results(search)
        filtered_results = []
        # select new publishment
        for r in tqdm(query_results):
            if r.updated.date() == self.date:
                paper = {
                        "Title": r.title,
                        "Summary": r.summary,
                        "Authors":[str(x) for x in r.authors],
                        "Link":str(r.links[0]),
                        "Source":"Arxiv"
                    }
                filtered_results.append(paper)
            else:
                break
        print(f"{len(filtered_results)} studies found in Arxiv.")
        self.source_results['Arxiv'] = filtered_results
        self.All_result.extend(filtered_results)
        
    def search_bioRxiv(self):
        keywords_list = self.configs['Keywords']
        target_date = self.date
        
        query = " OR ".join(keywords_list)
        query = query.replace(" ","%252B")

        base_url = "https://www.biorxiv.org"
        search_url = f"{base_url}/search/{query}%20numresults%3A{self.configs['max_num_per_source']}%20sort%3Apublication-date%20direction%3Adescending"
        print(f"query url as:\n {search_url}")

        try:
            # 发送 HTTP 请求
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            response = requests.get(search_url,headers)
            response.raise_for_status()

            # 使用 BeautifulSoup 解析 HTML
            soup = BeautifulSoup(response.text, "html.parser")

            # 提取文章链接
            links = []
            results = soup.find_all("li", class_="search-result", limit=self.configs['max_num_per_source'])
            for result in results:
                doi_span = result.find("span", class_='highwire-cite-metadata-doi')
                if doi_span:
                    doi_text = doi_span.get_text(strip=True)
                    link = doi_text.split("doi:")[-1].strip()  # 提取 "doi:" 后的部分
                else:
                    link = "No Link"
                links.append(link)
        except Exception as e:
            print(f"Error while scraping BioRxiv: {e}")
            return []

        # 检查每一篇文章的更新时间，提取信息
        filtered_results = []
        for link in tqdm(links):
            time.sleep(1)
            try:
                response = requests.get(link.replace("https://doi.org/","https://api.biorxiv.org/details/biorxiv/"),headers)
                response.raise_for_status()
                response = eval(response.text)

                # Extract date of update/publication
                update_date = datetime.strptime(response['collection'][0]['date'], "%Y-%m-%d").date()

                if update_date == target_date:
                    paper = {
                            "Title": response['collection'][0]['title'],
                            "Summary": response['collection'][0]['abstract'],
                            "Authors":response['collection'][0]["authors"].split('; '),
                            'Author_corresponding': response['collection'][0]["author_corresponding"],
                            'Author_corresponding_institution': response['collection'][0]["author_corresponding_institution"],
                            "Link":link,
                            "Source":"bioRxiv"
                        }
                    filtered_results.append(paper)
                elif update_date < target_date:
                    break
            except Exception as e:
                print(f"{link}:\n Error: {e}")

        print(f"{len(filtered_results)} studies found in BioRxiv.")
        self.source_results['bioRxiv'] = filtered_results
        self.All_result.extend(filtered_results)
        
        
    def search_PubMed(self):
        keywords_list = self.configs['Keywords']
        target_date = self.date
        
        query = " OR ".join(keywords_list)
        print(f"query condition as:\n {query}")

        query = "(" + query
        query += ")"    
        query += f'AND (\"{self.date.strftime("%Y/%m/%d")}\"[Date - Publication] : \"{(datetime.today().date()+ timedelta(days = 1)).strftime("%Y/%m/%d")}\"[Date - Publication])'

        try:
            fetch = PubMedFetcher()
            Entrez.email = ''
            pmids = Entrez.read(Entrez.esearch(db='pubmed', 
                                    sort="pub date", 
                                    retmax=self.configs['max_num_per_source'],
                                    retmode='xml', 
                                    term=query))['IdList']
            filtered_results = []
            for pmid in tqdm(pmids):
                info = fetch.article_by_pmid(pmid)
                paper = {
                        "Title": info.title,
                        "Summary": info.abstract,
                        "Authors":info.authors,
                        "Link":info.url,
                        "Source":info.journal
                        }
                filtered_results.append(paper)

            print(f"{len(filtered_results)} studies found in PubMed.")
            self.source_results['PubMed'] = filtered_results
            self.All_result.extend(filtered_results)

        except Exception as e:
            print(f"Error while scraping PubMed: {e}")
    ##############################
    
    
    ########## functions about query resource ##########
    def _remove_duplicates(self):
        """把不同来源查到的同样的文献删去"""
        dict_list = self.All_result
        seen_titles = set()
        result = []
        for item in dict_list:
            title = item.get("Title")
            if title not in seen_titles:
                seen_titles.add(title)
                result.append(item)
        print(f"# Paper after removing duplicates: {len(result)}")
        self.All_result = result
        
    def _categorize_by_title(self):
        """根据题目构建dict，便于后面的分类"""
        categorized_data = {}
        for item in self.All_result:
            title = item.get("Title")
            if title:
                if title not in categorized_data:
                    categorized_data[title] = []
                categorized_data[title].append(item)
        self.All_result_dict =  categorized_data
    
    
    def _cluster_by_topic(self):
        """用大预言模型（通义千问）根据文章的标题（加摘要太贵了）对文章分类，并总结出话题"""
        try:
            client = OpenAI(
                api_key=self.configs['DASHSCOPE_API_KEY'], 
                base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
            )
            # 提出需求
            completion = client.chat.completions.create(
            model="qwen-plus", 
            messages=[
                {'role': 'system', 'content': 'You are a helpful assistant to organize scientific papers.'},
                {'role': 'user', 'content': 
                 f"""Here are a list of paper titles, clustering these paper with their titles and give each cluster a suitable topic. Control the number of topics and avoid to put one paper in each topic. Only return me a structral dict like:
                 dict("Topic1":["Title 1", "Title 2", ...],
                 "Topic2":["Title 1", "Title 2", ...],
                 ...).Be Careful with the punctuation marks in the original title. Make sure the titles in the output is exactly same with the input. \nHere is the list {list(self.All_result_dict.keys())}"""}
            ])
            # organize the output of the LLM into a dict
            topic_dict_string = completion.to_dict()['choices'][0]['message']['content']
            topic_dict_string = topic_dict_string.strip("```")
            topic_dict_string = topic_dict_string.strip("python")
            topic_dict = eval(topic_dict_string)
            print(f"Clustered {len(self.All_result_dict)} Papers into {len(topic_dict)} Topics ")
            self.topic_dict = topic_dict

        except Exception as e:
            print(f"Error while using DASHSCOPE: {e}")
            
            
    def Process_collected_papers(self):
        self._remove_duplicates()
        self._categorize_by_title()
        self._cluster_by_topic()
            
    ##############################
    
    
    ########## functions about News ##########
        
    def _get_url_by_date(self, base_url, target_date):
        """
        获取特定日期对应的网址
        :param base_url: 目标网站的基础URL
        :param target_date: 目标日期 (格式：Dec 5, 2024)
        :return: 对应日期的网址或提示信息
        """
        try:
            # 请求网页内容
            response = requests.get(base_url)
            response.raise_for_status()  # 确保请求成功

            # 解析网页
            soup = BeautifulSoup(response.text, 'html.parser')

            # 找到包含目标日期的 div，并获取其父 <a> 标签的 href
            date_divs = soup.find_all('div', text=target_date)  # 匹配日期的 <div> 标签

            for div in date_divs:
                parent_a = div.find_parent('a')  # 获取父级 <a> 标签
                if parent_a and 'href' in parent_a.attrs:
                    return "https://www.aibase.com" + parent_a['href']  # 拼接完整链接

            return None
        except requests.RequestException as e:
            return None
        except Exception as e:
            return None
        
    def _get_first_ai_daily_article_from_news(self, base_url):
        """
        在 /news 页面获取第一篇标题包含 'AI Daily' 的文章链接
        :param base_url: 新闻页面的基础URL
        :return: 第一篇标题包含 'AI Daily' 的文章链接
        """
        try:
            # 请求网页内容
            response = requests.get(base_url)
            response.raise_for_status()  # 确保请求成功

            # 解析网页内容
            soup = BeautifulSoup(response.text, 'html.parser')

            # 查找标题包含 'AI Daily' 的 <h3> 标签
            articles = soup.find_all('h3', text=lambda x: x and "AI Daily" in x)

            if articles:
                first_article = articles[0]  # 获取第一篇文章的 <h3>
                parent_a = first_article.find_parent('a')  # 获取父级 <a> 标签
                if parent_a and 'href' in parent_a.attrs:
                    return "https://www.aibase.com" + parent_a['href']  # 拼接完整链接

            return None
        except requests.RequestException as e:
            return None
        except Exception as e:
            return None
        
    def _extract_titles_from_page(self, url):
        """
        从网页中提取所有 <strong> 元素的内容
        :param url: 目标网页URL
        :return: 标题列表
        """
        try:
            # 请求网页内容
            response = requests.get(url)
            response.raise_for_status()  # 确保请求成功

            # 解析网页内容
            soup = BeautifulSoup(response.text, 'html.parser')

            # 查找所有符合条件的 <strong> 元素
            strong_elements = soup.find_all('strong', class_="")

            # 提取每个 <strong> 标签中的文本
            titles = [element.get_text(strip=True) for element in strong_elements]
            cleaned_titles = [title for title in titles if title and title != "Click to Learn More"]
            return cleaned_titles
        except requests.RequestException as e:
            print(f"请求失败: {e}")
            return []
        except Exception as e:
            print(f"解析失败: {e}")
            return []
    
    def get_News(self):
        news_url = self._get_url_by_date("https://www.aibase.com/daily", self.date.strftime('%b %d, %Y'))
        if news_url is None:
            news_url = self._get_first_ai_daily_article_from_news("https://www.aibase.com/news")
        
        self.News_link = news_url
        if news_url is None:
            self.News = []
            print("News not available today.")
        else:
            self.News = self._extract_titles_from_page(news_url) 
    ##############################
    
    ########## functions about construct report ##########
    def construct_output(self):
        md_content = f"""
<!-- Logo 图片和日期 -->
<p align="center" style="display: flex; align-items: center; justify-content: center;">
  <img src="./ScholarDaily.png" alt="Logo" width="200" style="margin-right: 20px;">
  <span style="font-size: 16px; color: gray;">{self.date.strftime("%Y/%m/%d")}</span>
</p>
"""
        # add papers
        md_content+="""
# Papers
---
"""
        for topic in self.topic_dict.keys():
            md_content += f"## {topic}\n"
            for idx, title in enumerate(self.topic_dict[topic]):
                try:
                    md_content += f"### {idx+1}. [{title}]({self.All_result_dict.get(title)[0]["Link"]})\n"
                    md_content += (f"* **Journal**: {self.All_result_dict.get(title)[0].get("Source", "")}\n")
                    md_content += ("* **Authors**:" + ", ".join(self.All_result_dict.get(title)[0].get("Authors", ""))+"\n")
                    if self.All_result_dict.get(title)[0].get("Summary", "") is None:
                        md_content += "Summary not available."
                    else:
                        md_content += self.All_result_dict.get(title)[0].get("Summary", "")
                    md_content += "\n\n"
                except Exception as e:
                    print(f"In paper {title}:\n found error{e}")

        # add news
        md_content+="""# AI News (excerpted from AIbase)
---
"""
        for idx,news in enumerate(self.News):
            md_content+=f"### {news}\n"

        md_content+=f"see more details on [AIbase]({self.News_link})" 

        md_content += "\n\n"



        ## end of the document
        md_content+="""---
<p align="center" style="font-family: \'Great Vibes\', cursive; font-size: 20px;">————————   Stay ahead, stay inspired.   ————————</p>"""

        md_content+="""
<!-- 固定在右下角的联系方式 -->
<div style="position: relative; margin-top: 30px; font-size: 14px; color: gray; text-align: center;">
    Contact：<br>
    📧 Email: ysionchen@gmail.com<br>
    🌐 github: <a href="https://github.com/cyxss/ScholarDaily" target="_blank" style="color: gray;">github.com/cyxss/ScholarDaily</a>
</div>"""
        
        
        html_content = markdown.markdown(md_content, extensions=['extra', 'toc', 'tables'])

        full_html = f"""

<!DOCTYPE html> <html lang="en"> <head> <meta charset="UTF-8"> <meta name="viewport" content="width=device-width, initial-scale=1.0"> <title>Markdown to HTML</title> <style> body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }} img {{ display: block; margin: 0 auto; }} p {{ text-align: center; }} </style> </head> <body> {html_content} </body> </html> """
        
        
        html_file_path = f"ScholarDaily_{self.date}.html" 
        with open(html_file_path, "w", encoding="utf-8") as file: 
            file.write(full_html)
    ##############################
    
    def produce(self):
        # searching in each source 
        for source in self.configs["Sources"]:
            if source in self.fuction_mapping.keys():
                self.fuction_mapping.get(source)()
            else:
                print(f"""{source} tool not define!
Please define the function and add it into the ScholarDaily.fuction_mapping as 'source_name': self.function_name
                      """)
                return -1
            
        # process the collected papers
        self.Process_collected_papers()
        
        # collect news
        self.get_News()
        
        # construct output
        self.construct_output()
        
        # logging
        
        with open(f"{self.configs.get("log_path","./")}/{self.date}.pkl", 'wb') as file:
            pickle.dump(self, file)
