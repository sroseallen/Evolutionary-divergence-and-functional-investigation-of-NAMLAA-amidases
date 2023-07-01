import urllib.request
import urllib.parse

parameters = {
    'seqdb':'uniprotrefprot',
    'seq':'VVLDPGHGGIDGGARGVTGILEKDVTLAFARALRDELQKGSHTIVALTRDSDIFLRLSERVKKAQEFDADLFISIHADTIDVHSLRGATVYTISDEASDAIAKSLAESENKVDLLDGLPKEESLELTDILLDLTRRETHAFSINFANNVVSNLSKSHINLINNPHRYADFQVLKAPDVPSVLIEIGYLSNKEDEKLLNNPQWRKQMAASIAYSIRQ'
}
enc_params = urllib.parse.urlencode(parameters)
enc_params = enc_params.encode('ascii')

#post the search request to the server
request = urllib.request.Request('https://www.ebi.ac.uk/Tools/hmmer/search/phmmer/', enc_params)
#results_url = urllib.request.urlopen(request)

# send GET request
#res_params = {
#    'output': 'json',
#    'range': '1,10'
#}

# add the parameters to your request for the results
#enc_res_params = urllib.parse.urlencode(res_params)
#enc_res_params = enc_res_params.encode('ascii')
#modified_res_url = results_url + enc_res_params

#results_request = urllib.request.Request(modified_res_url)
#data = urllib.request.urlopen(results_request)

with urllib.request.urlopen(request) as f:

    print(f.read().decode('utf-8'))