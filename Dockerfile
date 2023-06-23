FROM python:3.10

WORKDIR /app

COPY requirements.txt .

RUN apt-get update

RUN pip install --no-cache-dir -r requirements.txt

COPY src ./src

CMD ["python", "src/main.py"]