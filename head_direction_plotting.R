#install.packages('RMySQL') #install the package if you don't have it already
library('ggplot2')
library('RMySQL') #load the necessary library

database = dbConnect(MySQL(), user='root', password='q88ul6014!', dbname='data', host='localhost') #connect to server, replace fields with your values

dbListTables(database) #show the tables in the database

rs = dbSendQuery(database, "SELECT Geometry,Genotype,Sessions,Tetrode,Unit,ID,MVL,peak_direction,angular_stability_binned, directional_information,avg_rate_spd_flt FROM all_cells WHERE Matching_squish = 1 AND Geometry = 'Open Field' AND MVL >= 0.2 ORDER BY Genotype
") #run a query as "rs". This particular query will grab all the data from the table `table`

open_field_data = fetch(rs, n=-1) #actually grab the data from the query "rs" and store it in the dataframe "data"

rt = dbSendQuery(database, "SELECT Geometry,Genotype,Sessions,Tetrode,Unit,ID,MVL,peak_direction,angular_stability_binned, directional_information,avg_rate_spd_flt FROM all_cells WHERE Matching_squish = 1 AND Geometry = 'Squish' AND MVL >= 0.2 ORDER BY Genotype")

squish_data = fetch(rt, n=-1)

angles = (open_field_data$peak_direction)
peak_direction = as.vector(open_field_data$peak_direction)

ggplot(aes(x = peak_direction, y = MVL), data = open_field_data) + geom_point()
#+ coord_polar(peak_direction)
  