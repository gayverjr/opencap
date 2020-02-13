#include <iostream>
#include <vector>
#include "overlap.h"
#include "BasisSet.h"
#include <armadillo>
#include "utils.h"
#include "transforms.h"
using namespace std;

int main(int argc, char **argv)
{
	BasisSet bs("N2.xyz","bas.bas");
	//code block for generating overlap matrix, re-use general principles elsewhere
	arma::mat Smat(bs.num_carts(),bs.num_carts());
	Smat.zeros();
	unsigned int row_idx = 0;
	for(auto&shell1:bs.basis)
	{
		unsigned int col_idx = 0;
		for(auto&shell2: bs.basis)
		{
			auto sub_mat = Smat.submat(row_idx,col_idx,
					row_idx+shell1.num_carts()-1,col_idx+shell2.num_carts()-1);
            shell_overlap(shell1,shell2,sub_mat);
            col_idx += shell2.num_carts();
		}
		row_idx += shell1.num_carts();
	}
	std::cout << "COEFF:" << get_coeff(4,0,2,0,2) << std::endl;
	//Smat.print();
	uniform_cart_norm(Smat,bs);
	//Smat.print();
	//cart2spherical(Smat,bs);
	//Smat.print();
	return 0;
}

