import matplotlib.dates as mdates
import matplotlib.ticker as ticker

def format_line_ticks(ax, year_range):
    year_range = abs(year_range)
    print "year_range = ", year_range
    if year_range < 0.0005:       # a few hours
        ax.xaxis.set_major_locator(mdates.HourLocator())
        ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=4))

        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%M"))
        print "ft: 0; %s" % year_range

    elif year_range < 1.1/365:       #about a day
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=5))
        ax.xaxis.set_minor_locator(mdates.HourLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter(""))
        print "ft: 0.5; years: %s" % year_range
        print "         days:  %s" % (year_range*365)

    elif year_range < 3.1/365:       #about three days
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
        ax.xaxis.set_minor_locator(mdates.HourLocator(interval=4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter(""))
        print "ft: 0.75; year: %s" % year_range
        print "          days: %s" % (year_range*365)

    elif year_range < 7.1/365:       #about a week
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter(""))
        print "ft: 0.85; year: %s" % year_range
        print "          days: %s" % (year_range*365)

    elif year_range < 0.05:       #about half a month

        ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter(""))
        print "ft: 1; %s" % year_range
        print "      days: %s" % (year_range*365)
    
    elif year_range < 0.18:       #about two months
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1, interval=2))
        ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=1))

        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        
        print "ft: 2; yrs: %s" % year_range
        print "      days: %s" % year_range*365.0
    
    
    elif year_range < 0.3:       #about 3.5 months
        ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1))
        ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=1))

        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        
        print "ft: 3; %s" % year_range
        
         

    elif year_range < 0.5:       #about 6 months
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))

        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_formatter(mdates.DateFormatter('%b %Y'))

        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
            
        print "ft: 4; %s" % year_range
    
    elif year_range < 1.1:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))

        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))

        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
        
        print "ft: 5; %s" % year_range

    elif year_range < 18.3:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_minor_locator(mdates.YearLocator(month=7))

        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("'%y"))

        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment('center')
        
        print "ft: 6; %s" % year_range
            
    else:
        ax.xaxis.set_major_locator(mdates.YearLocator(10))
        ax.xaxis.set_minor_locator(mdates.YearLocator(10))

        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    
        print "ft: 7; %s" % year_range
