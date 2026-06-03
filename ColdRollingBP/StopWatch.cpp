//
// Created by wshikui on 2023/10/17.
//
#include "StopWatch.h"

Stopwatch::Stopwatch(const char* name) :
	m_name(name)
{
	start();
}

Stopwatch::Stopwatch(const std::string& name) :
	m_name(name.c_str())
{
	start();
}


Stopwatch::~Stopwatch()
{
	if (!m_stopped)
	{
		stop();
	}
}


double Stopwatch::elapsed() const
{
	auto end = std::chrono::high_resolution_clock::now();

	return std::chrono::duration_cast<std::chrono::duration<double>>(
		end - m_start).count();
}


double Stopwatch::lap()
{
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(
		end - m_last_time).count();
	m_last_time = end;
	return duration;
}


void Stopwatch::start()
{
	m_stopped = false;
	m_start = std::chrono::high_resolution_clock::now();
	m_last_time = m_start;
}


double Stopwatch::stop()
{
	auto end = std::chrono::high_resolution_clock::now();

	m_stopped = true;

	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(
		end - m_start).count();

	return duration;
}