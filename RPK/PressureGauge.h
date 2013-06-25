/**
 * PressureGauge.h - pressure guage class declaration
 * Written by Alexanfer Bass, modified by Max I. Fomitchev-Zamilov and
 * Spenser Bauman
 */

#ifndef __RPK_PRESSURE_GAUGE_H__
#define __RPK_PRESSURE_GAUGE_H__

#include <boost/array.hpp>
#include <boost/circular_buffer.hpp>
#include <deque>

const size_t PGAUGE_ORDER = 2;

typedef std::pair<double, double> TimeSeries;

class PressureGauge
{
public:

#ifdef WINDOWED_PGAUGE
    typedef std::deque<TimeSeries> BufferType;
#else
    typedef boost::circular_buffer<TimeSeries> BufferType;
#endif

    /**
     * Default Constructor. Creates a zeroed pressure gauge.
     */
	PressureGauge();

    /**
     * Add a sample to the pressure gauge.
     *
     * @param pressure A sample of the pressure to add to the gauge.
     */
	void AddSample(const double time, const double pressure);


    /**
     * Computes the time of the pressure sample that will be produced by
     * a call to the <code>SampleGauge</code> function. This is needed before
     * the call, as the time is used to determine the appropriate radius to
     * call <code>SampleGauge</code> with.
     */
    double PressureTime() const;

    /**
     * Get the pressure at the given time. This given time value is assumed
     * to be greater than the previous call.
     *
     * @param currentTime The current simulation time in simulation units.
     * @param pressure The location at which the computed pressure is stored.
     * @param dpdt The location at which the computed pressure derivative is
     *             stored.
     *
     * @return The time at which the sample is said to occur.
     */
	double SampleGauge(const double radius);

    /**
     * Compute the pressure at the given time. This makes use of the
     * coefficients of the polynomial calculated by the regression.
     */
    double ComputePressure(const double currentTime,
                           const double radius,
                           const double drdt) const;

    /**
     * Computes the rate of change of the pressure using the coefficients
     * of the polynomial calculated from the regression.
     */
    double ComputeDpdt(const double currentTime,
                       const double radius,
                       const double dpdt) const;

    /**
     * Serialize the pressure gauge to disk.
     *
     * @param file The file to save the pressure gauge to.
     * @param save Flag to control saving vs loading. True saves the file to
     *             disk while false causes the Serialize function to read from
     *             the given file.
     */
	void Serialize(FILE* file, bool save = true);

    template <typename Writer>
    void Serialize(Writer& pt) const;

    // void Read(const boost::property_tree::ptree& tree);

    bool IsReady() const;

    bool CanSample() const;

    void SetPressureStats(const double t, const double p, const double dpdt);

private:

    // Circular buffer of pressure samples.
    BufferType Buffer;

    boost::array<double, PGAUGE_ORDER> Coefficients;

    double SampleRadius;
    double SampleTime;

    double LastPressure;
    double LastDpdt;
    double LastTime;
    double CurrentPressure;
    double CurrentDpdt;
    double CurrentTime;

    bool Ready;
};

template <typename Writer>
void PressureGauge::Serialize(Writer& writer) const
{
    writer.StartObject();

    writer.String("LastPressure")    , writer.Double(LastPressure);
    writer.String("LastDpdt")        , writer.Double(LastDpdt);
    writer.String("LastTime")        , writer.Double(LastTime);
    writer.String("CurrentPressure") , writer.Double(CurrentPressure);
    writer.String("CurrentDpdt")     , writer.Double(CurrentDpdt);
    writer.String("CurrentTime")     , writer.Double(CurrentTime);

    // Write the buffer as an array of pairs. Each element of the array is
    // an object with two values: "first" and "second".
    writer.String("Buffer");
    writer.BeginArray();
    std::for_each(Buffer.begin(), Buffer.end(), [&writer](const TimeSeries& ts)
            {
                writer.BeginObject();
                writer.String("first"),  writer.Double(ts.first);
                writer.String("second"), writer.Double(ts.second);
                writer.EndObject();
            });

    writer.EndArray();

    // End the pressure gauge object.
    writer.EndObject();
}

#endif
