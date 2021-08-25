#include "include/user_interface.h"


std::random_device rd;
std::mt19937::result_type seed = rd() ^ (
	(std::mt19937::result_type)
	std::chrono::duration_cast<std::chrono::seconds>(
		std::chrono::system_clock::now().time_since_epoch()
		).count() +
	(std::mt19937::result_type)
	std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::high_resolution_clock::now().time_since_epoch()
		).count());
std::mt19937_64 gen(seed);


int main(const int argc, char* argv[]) {
	std::unique_ptr<user_interface> intface = std::make_unique<isingUI::ui>(argc, argv);
	intface->make_sim();
	return 0;
}