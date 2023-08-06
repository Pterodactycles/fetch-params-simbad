from querysimbad import querysimbad as qs

### query a user entered list of stars and plot the altitude for a given night

# query = qs()
# query.plot_alt(["hd_147584","hd_93486"], '2023-08-06')



### query the same as above but for one star

query = qs()
query.plot_alt("hd_147584", '2023-08-06')



### query the same but plot for the whole year

# query = qs()
# query.plot_alt_year(["hd_147584","hd_93486"])



### lets do this with 1000 stars, see that it doesn't break
# inputfilepath = "starlist.txt"
# with open(inputfilepath, newline="\n") as f:
#     inputfiledata = f.read().splitlines()

# query = qs()
# query.plot_alt_year(inputfiledata, legend=False)



### create a file of parameters from simbad for a list of stars

# query = qs()
# query.query_star_list("starlist.txt", batchsize=300)



### retrieve the parameters from simbad for a single star, or a couple of stars

# query = qs()
# onestar = query.retrieve_data("hd_147584")
# twostars = query.retrieve_data(["hd_147584","hd_93486"])


