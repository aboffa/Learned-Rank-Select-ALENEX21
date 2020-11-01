#include <iostream>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "index_reader.hpp"
#include "datasets_generation_utils.hpp"

template<typename TypeIn>
void write_data_binary(std::vector<TypeIn> *vec, std::string &filename, bool first_is_size = true,
                       size_t max_size = std::numeric_limits<size_t>::max() - 1) {
    try {
        auto n = vec->size();
        std::ofstream fout(filename, std::ios::out | std::ios::binary);
        if (first_is_size)
            fout.write((char *) &(n), sizeof(size_t));
        fout.write((char *) vec->data(), n * sizeof(TypeIn));
        fout.close();
    }
    catch (std::ios_base::failure &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << std::strerror(errno) << std::endl;
        exit(1);
    }
}

int main(int argc, char *argv[]) {
    std::cout << std::scientific;
    if (argc < 2) {
        std::cout << "USAGE: Generate_datasets <operation> <file_1> <file_2> ... <file_n> \n"
                     "Operations: \n"

                     "  - write_integers_from_BWT <- files must be Burrow Wheeler Transform file (*.bwt) obtained using"
                     " bwtdisk library (https://people.unipmn.it/manzini/bwtdisk/).\n"

                     "  - write_integers_from_genome <- files must contain just 'A', 'C', 'G' and 'T' characters.\n"

                     "  - write_integers_from_GOV2 <- files must be Inverted Indexes (for example GOV2) that are"
                     " readable using the class IndexReader in the file index_reader.hpp.\n"

                     "  - write_google_ngram_v1 <- files must be a Google 5gram file like: googlebooks-eng-all-5gram-20090715-xxx.csv"

                     << std::endl;

        return -1;
    }

    std::string path_to_save = "datasets/";

    if (strcmp(argv[1], "write_integers_from_BWT") == 0) {
        using type = uint32_t;
        for (int i = 2; i < argc; ++i) {
            std::vector<WT_char_positions> characteristic_vectors(256);
            std::string file = std::string(argv[i]);
            std::ifstream myfile(file, std::ios_base::binary | std::ios_base::in | std::ios_base::ate);

            int64_t num_bytes = 0;
            if (myfile.is_open()) {
                num_bytes = myfile.tellg();
                myfile.seekg(17);

                for (auto j = 0; j < 256; ++j) {
                    characteristic_vectors[j] = WT_char_positions(new std::vector<uint32_t>(), (char) j);
                    characteristic_vectors[j].vec->reserve(num_bytes);
                }

                for (auto num_char = 17; num_char < num_bytes - 17; ++num_char) {
                    char c = 0;
                    myfile.read(&c, 1);
                    if (isalnum(c) || isspace(c) || ispunct(c)) { // ecluding control character
                        characteristic_vectors[(size_t) (c)].vec->push_back(num_char);
                    }
                }
            }
            std::sort(characteristic_vectors.begin(), characteristic_vectors.end(),
                      [](const WT_char_positions &a, const WT_char_positions &b) {
                          return a.vec->size() > b.vec->size();
                      });

            for (auto j = 0; j < 2; ++j) {
                std::string filename = path_to_save + std::to_string(j) + "_WT_" + file;
                write_data_binary<type>(characteristic_vectors[j].vec, filename, true);
            }
            float ratio = static_cast<float>(characteristic_vectors[2].vec->size()) / num_bytes;
            float wanted_ratio = ratio;
            int times = 2;
            size_t index = 2;
            while (times < 8 && index < characteristic_vectors.size()) {
                ratio = static_cast<float>(characteristic_vectors[index].vec->size()) / num_bytes;
                if (ratio < wanted_ratio * 1.1) {
                    wanted_ratio /= 2;
                    std::string filename = path_to_save + std::to_string(times) + "_WT_" + file;
                    write_data_binary<type>(characteristic_vectors[index].vec, filename, true);
                    times++;
                }
                index++;
            }
        }
    }

    if (strcmp(argv[1], "write_integers_from_genome") == 0) {
        using type = uint32_t;
        srand(21);
        for (int i = 2; i < argc; ++i) {
            std::vector<WT_char_positions> characteristic_vectors(256);
            std::string file = std::string(argv[i]);
            std::ifstream myfile(file, std::ios_base::binary | std::ios_base::in | std::ios_base::ate);

            int64_t num_bytes = 0;
            if (myfile.is_open()) {
                num_bytes = myfile.tellg();
                myfile.seekg(0);

                for (auto j = 0; j < 256; ++j) {
                    characteristic_vectors[j] = WT_char_positions(new std::vector<uint32_t>(), (char) j);
                    characteristic_vectors[j].vec->reserve(num_bytes);
                }

                for (auto num_char = 0; num_char < num_bytes; ++num_char) {
                    char c = 0;
                    myfile.read(&c, 1);
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                        std::cout << "Wrong file format" << std::endl;
                        exit(1);
                    }
                    characteristic_vectors[(size_t) (c)].vec->push_back(num_char);
                }

                for (auto j = 0; j < 256; ++j) {
                    characteristic_vectors[j].vec->shrink_to_fit();
                }
            }
            std::sort(characteristic_vectors.begin(), characteristic_vectors.end(),
                      [](const WT_char_positions &a, const WT_char_positions &b) {
                          return a.vec->size() > b.vec->size();
                      });

            for (auto j = 0; j < 4; ++j) {
                std::string filename = path_to_save + std::to_string(j) + "_" + characteristic_vectors[j].c + "_" + file;
                write_data_binary<type>(characteristic_vectors[j].vec, filename, true);
            }

            for(auto index = 0; index < 4; index++) {
                size_t n = 0;
                double ratio = 0;
                size_t u = num_bytes;
                std::vector<type> old_data;
                old_data.swap(*(characteristic_vectors[index].vec));
                std::vector<type> new_data;

                for (int h = 0; h < 2; ++h) {
                    new_data.clear();
                    new_data.reserve(old_data.size());
                    for (unsigned int j : old_data) {
                        if (rand() % 5 == 0) {
                            new_data.push_back(j);
                        }
                    }
                    old_data.clear();
                    new_data.shrink_to_fit();
                    n = new_data.size();
                    ratio = static_cast<double>(n) / u;
                    std::cout << n << "," << u << "," << ratio << std::endl << std::flush;
                    std::string filename =
                            path_to_save + std::to_string(h + 4) + "_" + characteristic_vectors[index].c + "_" + file;
                    write_data_binary<type>(&new_data, filename, true);
                    old_data.swap(new_data);
                }
            }
        }
    }


    if (strcmp(argv[1], "write_integers_from_GOV2") == 0) {
        using type = uint32_t;
        std::string file = std::string(argv[2]);
        IndexReader ir(file);
        std::vector<std::vector<type>> invertedLists;
        std::vector<type> buffer;
        while (ir.load_integers(buffer)) {
            invertedLists.push_back(buffer);
            buffer.clear();
        }

        std::sort(invertedLists.begin(), invertedLists.end(),
                  [](const std::vector<type> &a, const std::vector<type> &b) {
                      return a.size() > b.size();
                  });

        for (auto i = 0; i < invertedLists.size(); ++i) {
            if (invertedLists[i].size() < 100000) {
                invertedLists.erase(invertedLists.begin() + i, invertedLists.end());
                break;
            }
        }

        auto times = 0;
        std::vector<std::string> directory_names;
        directory_names.emplace_back("100K-1M/");
        directory_names.emplace_back("1M-10M/");
        directory_names.emplace_back("10M-/");

        for (auto data = invertedLists.begin(); data != invertedLists.end(); ++data) {
            times++;
            std::cout << std::endl;
            std::string filename;
            switch (data->size()) {
                case 100000 ... 1000000:
                    filename = "datasets/GOV2/" + directory_names[0] + std::to_string(times) + "_" + file;
                    write_data_binary<type>(&(*data), filename, true);
                    break;
                case 1000001 ... 10000000:
                    filename = "datasets/GOV2/" + directory_names[1] + std::to_string(times) + "_" + file;
                    write_data_binary<type>(&(*data), filename, true);
                    break;
                default:
                    filename = "datasets/GOV2/" + directory_names[2] + std::to_string(times) + "_" + file;
                    write_data_binary<type>(&(*data), filename, true);
                    break;
            }

        }
    }

    if (strcmp(argv[1], "write_google_ngram_v1") == 0) {
        using type = uint32_t;
        std::vector<std::string> data;
        std::vector<std::string> filenames;
        filenames.reserve(argc);
        for (int i = 2; i < argc; ++i) {
            filenames.emplace_back(argv[i]);
        }
        std::string base_path;

        // template parameter is the presence of duplicates
        auto stats = load_data_from_Google_nGram_v1<false>(data, base_path, filenames);

        std::cout << "average LCP:  " << stats->averageLCP << std::endl;
        std::cout << "average length:  " << stats->averageLength << std::endl;
        std::cout << "max size:  " << stats->maxSize << std::endl;
        std::cout << "min size:  " << stats->minSize << std::endl;
        std::cout << "size:  " << stats->size << std::endl;
        std::cout << "nchars:  " << stats->chars.size() << std::endl;

        std::ofstream myfile;
        myfile.open("ALL_google-eng-all-5gram.txt");
        for (auto i = 0; i < data.size(); ++i) {
            myfile.write(data[i].c_str(), data[i].size());
            myfile.write("\n", 1);
        }
    }
}
