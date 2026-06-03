#pragma once
#include <chrono>
#include <string>

class Stopwatch
{
public:

	Stopwatch() = delete;
	/**
	 * @brief default constructor, copy constructor, move constructor, copy
	 * assignment, move assignment.
	*/
	Stopwatch(const Stopwatch& other) = default;
	Stopwatch(Stopwatch&& other) = default;
	Stopwatch& operator=(const Stopwatch& other) = default;
	Stopwatch& operator=(Stopwatch&& other) = default;

	/**
	 * @brief constructor.
	 * @param: const char*: scope name.
	*/
	// Stopwatch(const char* name, Fn &&func);
	explicit Stopwatch(const char* name);
	explicit Stopwatch(const std::string& name);

	/**
	 * @destructor
	*/
	~Stopwatch();

	/**
	 * @brief Elapsed time in seconds since the start of the stopwacth.
	 * @return double: time in seconds.
	*/
	[[nodiscard]] double elapsed() const;

	/**
	 * @brief Elapsed time in seconds since the last call of lap.
	 * @return double: time in seconds.
	*/
	double lap();

	/**
	 * @brief Starts (or re-starts) the timer.
	*/
	void start();

	/**
	 * @brief Stops the timer.
	*/
	double stop();

private:

	// using time_pt = std::chrono::time_point<std::chrono::system_clock>; old version
	using time_pt = std::chrono::time_point<std::chrono::high_resolution_clock>;

	/**
	 * @brief scope name.
	*/
	const char* m_name;

	/**
	 * @brief start time point.
	*/
	time_pt m_start;

	/**
	 * @brief.
	*/
	time_pt m_last_time;

	/**
	 * @brief.
	*/
	bool m_stopped{};

};
