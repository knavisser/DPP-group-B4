import json, pprint, requests, textwrap
host = 'http://localhost:8998'
data = {'kind': 'pyspark'}
headers = {'Content-Type': 'application/json'}
session_url = 'http://localhost:8998/sessions/0'
statements_url = session_url + '/statements'
r = requests.post(host + '/sessions', data=json.dumps(data), headers=headers)
r.json()

{u'id': 1, u'state': u'idle'}

data = {
  'code': textwrap.dedent("""
    import random
    NUM_SAMPLES = 100000
    def sample(p):
      x, y = random.random(), random.random()
      return 1 if x*x + y*y < 1 else 0

    count = sc.parallelize(xrange(0, NUM_SAMPLES)).map(sample).reduce(lambda a, b: a + b)
    print "Pi is roughly %f" % (4.0 * count / NUM_SAMPLES)
    """)
}

r = requests.post(statements_url, data=json.dumps(data), headers=headers)
pprint.pprint(r.json())

{u'id': 12,
u'output': {u'data': {u'text/plain': u'Pi is roughly 3.136000'},
            u'execution_count': 12,
            u'status': u'ok'},
u'state': u'running'}