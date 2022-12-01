#spark-submit \
#  --master spark://spark-master:7077 \
#  /opt/homebrew/Cellar/apache-spark/3.3.1/libexec/examples/src/main/python/pi.py \
#  2


spark-submit \
  --master local \
  /opt/homebrew/Cellar/apache-spark/3.3.1/libexec/examples/src/main/python/pi.py \
  100


