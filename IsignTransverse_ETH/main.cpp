#include "include/user_interface.h"

int main(const int argc, char* argv[]) {
	std::unique_ptr<user_interface> intface = std::make_unique<isingUI::ui>(argc, argv);
	intface->make_sim();
	return 0;
}