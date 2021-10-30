#pragma once
#include<thread>
#include<vector>
#include<queue>
#include <condition_variable>
#include<mutex>
#include<functional>
#include<memory>
#include<future>
#include<iostream>
#include <type_traits>
#include <numeric>   

//Thread safe queue from C++ Concurrency in Action by Anthony Williams
//https://www.manning.com/books/c-plus-plus-concurrency-in-action-second-edition

template<typename T>
class threadsafe_queue
{
private:
	mutable std::mutex mut;
	std::queue<std::shared_ptr<T> > data_queue;
	std::condition_variable data_cond;
public:
	threadsafe_queue()
	{}
	void wait_and_pop(T& value)
	{
		std::unique_lock<std::mutex> lk(mut);
		//https://www.youtube.com/watch?v=eh_9zUNmTig watch this to understand
		data_cond.wait(lk, [this] {return !data_queue.empty(); });
		value = std::move(*data_queue.front());
		data_queue.pop();
	}

	bool try_pop(T& value)
	{
		std::lock_guard<std::mutex> lk(mut);
		if (data_queue.empty())
			return false;
		value = std::move(*data_queue.front());//dereference the pointer
		data_queue.pop();
		return true;
	}

	std::shared_ptr<T> wait_and_pop()
	{
		std::unique_lock<std::mutex> lk(mut);
		data_cond.wait(lk, [this] {return !data_queue.empty(); });
		std::shared_ptr<T> res = data_queue.front();
		data_queue.pop();
		return res;
	}

	std::shared_ptr<T> try_pop()
	{
		std::lock_guard<std::mutex> lk(mut);
		if (data_queue.empty())
		{
			return std::shared_ptr<T>();
		}
		std::shared_ptr<T> res = data_queue.front();
		data_queue.pop();
		return res;
	}

	void push(T new_value)
	{
		std::shared_ptr<T> data(std::make_shared<T>(std::move(new_value)));
		std::lock_guard<std::mutex> lk(mut);
		data_queue.push(data);
		data_cond.notify_one();
	}

	bool empty() const
	{
		std::lock_guard<std::mutex> lk(mut);
		return data_queue.empty();
	}


};

//extends thread guard to vector of threads
class join_threads
{
	std::vector<std::thread>& _threads;

public:
	explicit join_threads(std::vector<std::thread>& threads)
		:_threads(threads)
	{}

	~join_threads()
	{
		for (auto& i : _threads)
		{
			if (i.joinable())
			{
				i.join();
			}
		}
	}
};


class function_wrapper
{
	struct impl_base {
		virtual void call() = 0;
		virtual ~impl_base() = default;
	};

	std::unique_ptr<impl_base>impl;

	template<typename F>
	struct impl_type : impl_base
	{
		F f;
		impl_type(F&& f_) : f(std::move(f_)) {}//constructor
		void call() { f(); }
	};

public:

	function_wrapper() = default;

	template<typename F>
	function_wrapper(F&& f) noexcept :
		impl(std::make_unique<impl_type<F>>(std::move(f)))
	{}


	void operator()() { impl->call(); }

	function_wrapper(function_wrapper&& other) :
		impl(std::move(other.impl))
	{
		other.impl = nullptr;
	}

	function_wrapper& operator=(function_wrapper&& other)
	{
		if (this != &other) {
			impl = std::move(other.impl);
			other.impl = nullptr;
		}

		return *this;
	}

	function_wrapper(const function_wrapper&) = delete;
	function_wrapper(function_wrapper&) = delete;
	function_wrapper& operator=(const function_wrapper&) = delete;
};


class thread_pool
{
	std::atomic_bool done;
	threadsafe_queue<function_wrapper> work_queue;
	std::vector<std::thread>threads;
	join_threads joiner;


	void worker_thread()
	{
		while (!done)
		{
			auto task = work_queue.try_pop();
			if (task != nullptr)
			{
				(*task)();
			}
			else
			{
				std::this_thread::yield();
			}
		}
	}

public:

	thread_pool() :
		done(false), joiner(threads)
	{
		unsigned const thread_count = std::thread::hardware_concurrency();

		try
		{
			for (unsigned int i = 0; i < thread_count; ++i)
			{
				threads.push_back(std::thread(&thread_pool::worker_thread, this));
			}
		}
		catch (...)
		{
			done = true;
			throw;
		}

	}

	~thread_pool()
	{
		//threads joining handled in constructor
		done = true;
	}

	/*template<typename FunctionType>
	auto submit(FunctionType f)->std::future<typename std::result_of<FunctionType()>>*/

	//attempt at perfect forwarding
	template<typename FunctionType, typename... Args>
	std::future<typename std::invoke_result_t<FunctionType, Args...>> submit(FunctionType&& f, Args&&... args)
	{
		using result_type = typename std::invoke_result_t<FunctionType, Args...>;
		/*std::packaged_task<result_type()> task(std::move(f));*/
		std::packaged_task<result_type()> task(std::bind(std::forward<FunctionType>(f), std::forward<Args>(args)...));
		std::future<result_type> res(task.get_future());
		work_queue.push(std::move(task));

		return res;
	}
};


//https://github.com/Toeplitz/threadpool example from here
int main()
{

	thread_pool tp;
	uint64_t N = 100;

	std::atomic_int ongoing_count(0);
	std::atomic_int completed_count(0);

	for (uint64_t n(0); n < N; ++n)
	{
		// Delay submitting to get clean looking output...
		std::this_thread::sleep_for(std::chrono::milliseconds(10));

		tp.submit([n, &ongoing_count, &completed_count]() {
			ongoing_count++;

			const uint32_t sleep_time = 1000;
			std::thread::id this_id = std::this_thread::get_id();
			std::cout << "[ " << n << "] sleeping for " << sleep_time << " in thread " << this_id << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time));
			std::cout << "[ " << n << "] exiting ..." << std::endl;

			completed_count++;
			ongoing_count--;
		});
	}

	while (1)
	{
		std::cout << "ongoing: " << ongoing_count << ", completed: " << completed_count
			<< " / " << N << std::endl;
		std::this_thread::sleep_for(std::chrono::milliseconds(1000));

		if (completed_count == N)
		{
			std::cout << "Finished: " << completed_count << std::endl;
			break;
		}
	}

}