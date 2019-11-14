for location in 'india' 'utrecht' 'stockholm' 'nairobi' 'manaus' 'kansas' 'chendu'
do
for env in 0
do
	for model in 15 # select the models you want to run 
	do
		for int in 1 #runs with or without interpolation
		do
			for parcels in 2000
			do
			echo $model $int $parcels './input/'$location$env'.nc .output/out'$location$env'_'$parcels'_'$int'_'$model'.nc'
		done
	done
done
done
done
