// Step 1: add header guards. What are header guards? Why are they necessary?
//So same code is not compiled multiple times. Pragma once can also be used

// Step 2: Create an empty template class MyVector. 

// Step 3: Implement the following methods:
// 1. void push_back(T data);
// 2. void insert(int index, T data);
// 3. T get (int index);
// 4. int size();
// 5. int capacity();
// 6. a destructor to release memory. 

// Self implementation of
// the Vector Class in C++

//#include <bits/stdc++.h>


#ifndef myvec
#define myvec


#include <iostream>

using namespace std;
template <typename T>



class vectorClass {
	

public:
T* arr;
	int capacity;
	int current;

	vectorClass()
	{
		arr = new T[1];
		capacity = 1;
		current = 0;
	}

	~ vectorClass()
	{
		delete [] arr;
	}

	void push_back(T data)
	{

		
		if (current == capacity) {
			T* temp = new T[capacity+1];

			for (int i = 0; i < capacity; i++) {
				temp[i] = arr[i];
			}

			delete[] arr;
			capacity++;
			arr = temp;
		}

		
		arr[current] = data;
		current++;
	}

	
	void insert(T data, int index)
	{

		if (index == capacity)
			push(data);
		else
			arr[index] = data;
	}

	T get(int index)
	{
		if (index < current)
			return arr[index];
	}


	int size() { return current; }

	int get_capacity() { return capacity; }

	void print()
	{
		for (int i = 0; i < current; i++) {
			cout << arr[i] << " ";
		}
		cout << endl;
	}
};


#endif