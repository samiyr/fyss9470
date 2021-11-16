#ifndef NUMBER_LIST_READER_H
#define NUMBER_LIST_READER_H

template<typename NumberType>
class NumberListReader {
public:
	NumberListReader(string input_file) {
		ifstream infile;
    	infile.open(input_file);

    	if (infile.is_open()) {
        	NumberType number;

        	while(infile >> number) {
            	numbers.push_back(number);
        	}
        }
	}

	const NumberType &operator()(typename std::vector<NumberType>::size_type index, typename std::vector<NumberType>::size_type offset = typename std::vector<NumberType>::size_type(0)) const {
		const auto requested_index = index + offset;
		auto actual_index = requested_index;
		const auto bounds = numbers.size();
		if (requested_index < 0) { actual_index = 0; }
		if (requested_index >= bounds) { actual_index = requested_index % bounds; }
		return numbers[actual_index];
	}
private:
	std::vector<NumberType> numbers;
};

#endif // NUMBER_LIST_READER_H