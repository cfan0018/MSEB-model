#!/usr/bin/env python
import cdsapi
 
c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'variable': ['geopotential', 'specific_humidity', 'temperature','u_component_of_wind', 'v_component_of_wind', 'vertical_velocity',],
        'product_type': 'monthly_averaged_reanalysis',
        'pressure_level': ['1000', '925', '850','700', '600', '500','400', '300', '250','200', '150', '100','70', '50', '30',],
        'year': ['1979','1980','1981','1982','1983','1984','1985','1986','1987','1988','1989','1990','1991','1992','1993','1994','1995','1996','1997','1998',
                 '1999','2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018',],
        'month': ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12',],
        'time': '00:00',
        'area': [30,-180,-30,180,],
        'grid': [3.75,3.75], 
    },  'era5prs.nc')

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': ['surface_latent_heat_flux', 'surface_net_solar_radiation', 'surface_net_thermal_radiation','surface_sensible_heat_flux', 'top_net_solar_radiation', 'top_net_thermal_radiation',],
        'year': ['1979','1980','1981','1982','1983','1984','1985','1986','1987','1988','1989','1990','1991','1992','1993','1994','1995','1996','1997','1998',
                 '1999','2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018',],
        'month': ['01', '02', '03','04', '05', '06','07', '08', '09','10', '11', '12',],
        'time': '00:00',
        'area': [30,-180,-30,180,],
        'grid': [3.75,3.75], 
    },  'era5sfc.nc')