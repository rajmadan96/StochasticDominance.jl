function ontosimplex!(probabilities::Vector{Float64})   # projects the probabilities onto the simplex
	for i= 1:length(probabilities)
		if probabilities[i] < 0.0 
			probabilities[i] = 0.0		# have them positive
		end
	end
	summ= sum(probabilities)
	if summ <= 0.0
		probabilities= fill!(probabilities, 1.0/ length(probabilities))
		summ= sum(probabilities)
		#@info "ontoSimplex: created probabilities." maxlog=1
	end
	if summ != 1
		probabilities./= summ
		tmpi= rand(1:length(probabilities))		# modify a random index
		summ= -sum(probabilities[1:tmpi-1])+ 1.0- sum(probabilities[tmpi+1:end])
		probabilities[tmpi]= max(0, summ)
		#@info "ontoSimplex: modified probabilities." maxlog=1
	end
	nothing		# modifies probabilities
end
