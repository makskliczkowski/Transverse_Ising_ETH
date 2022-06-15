
namespace thomas{
    
    inline
    auto tridiagonal_algorithm(
        arma::vec subdiagonal,      //<! one below main diagonal elements in matrix
        arma::vec diagonal,         //<! diagonal elements
        arma::vec superdiagonal,    //<! one above diagonal elements in matrix
        arma::vec result            //<! resulting vector
    ) -> arma::vec 
    {
        if(diagonal.size() != result.size()){
		    std::cout << "Incompatible dimensions: " << diagonal.size() << "vs.\t" << result.size() << std::endl;
		    assert(false);
	    }
        const size_t N = diagonal.size();
        arma::vec output(N, arma::fill::zeros);
        for(int i = 0; i < N -1; i++){
            const double w = subdiagonal(i) / diagonal(i);
            diagonal(i + 1) -= w * superdiagonal(i);
            result(i + 1) -= w * result(i);
        }
        output(N - 1) = result(N - 1) / diagonal(N - 1);
        for(int i = N - 2; i >=0; i--)
            output(i) = ( result(i) - superdiagonal(i) * output(i + 1) ) / diagonal(i);
        
        return output;
    }
}