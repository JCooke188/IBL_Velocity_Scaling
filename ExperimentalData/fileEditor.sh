# Created by Justin 
# Purpose is to take the files from the dune field calcs, remove the headers, and then change the file names


for file in x*comp*
do grep -v "#" "$file" > temp && mv temp "$file"
done

for file in *README
do grep -v "#" "$file" > temp && mv temp "$file"
done

