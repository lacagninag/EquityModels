﻿using System;
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using Fairmat.Math;
using Fairmat.Statistics;

namespace Dupire
{
    public partial class DupireEstimator : IEstimatorEx, IIntegrable
    {
        private EstimationResult QuantLibEstimate(InterestRateMarketData Mdataset, CallPriceMarketData Hdataset)
        {
            EquityCalibrationData HCalData = new EquityCalibrationData(Hdataset, Mdataset);

            bool hasArbitrage = HCalData.HasArbitrageOpportunity();
            if (hasArbitrage)
                Console.WriteLine("Market data contain arbitrage opportunity");

            this.r = new DVPLDOM.PFunction(null);
            this.q = new DVPLDOM.PFunction(null);
            this.r.Expr = (double[,])ArrayHelper.Concat(Mdataset.ZRMarketDates.ToArray(), Mdataset.ZRMarket.ToArray());
            this.q.Expr = (double[,])ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.DividendYield.ToArray());
            this.r.Parse(null);
            this.q.Parse(null);

            for (int i = 0; i < Hdataset.Volatility.R; i++)
            {
                for (int j = 0; j < Hdataset.Volatility.C; j++)
                {
                    Hdataset.Volatility[i, j] = Hdataset.Volatility[i, j] * Hdataset.Volatility[i, j] * Hdataset.Maturity[i];
                }
            }

            IFunction impVol = FitImplVolModel(Hdataset);

            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);
            prj.Symbols.Add(impVol);
            // doc.WriteToXMLFile("impVol.fair");

            // todo: spostare nei settings
            int nmat = 100;
            int nstrike = 100;
            double lastMat = Hdataset.Maturity[Hdataset.Maturity.Length - 1];
            double lastStr = Hdataset.Strike[Hdataset.Strike.Length - 1];
            Vector locVolMat = Vector.Linspace(0.0, lastMat, nmat);
            Vector locVolStr = Vector.Linspace(0.0, lastStr, nstrike);
            Matrix locVolMatrix = new Matrix(nmat, nstrike);
            double t, dt, forwardValue, y, dy, strike, strikep, strikem, w, wp, wm, dwdy;
            double d2wdy2, den1, den2, den3, strikept, strikemt, wpt, wmt, dwdt;
            Integrate integrate = new Integrate(this);

            for (int i = 0; i < nmat; i++)
            {
                t = locVolMat[i];
                forwardValue = Hdataset.S0 * Math.Exp(integrate.AdaptLobatto(0.0, t));
                for (int j = 0; j < nstrike; j++)
                {
                    strike = locVolStr[j];
                    y = Math.Log(strike / forwardValue);
                    dy = ((Math.Abs(y) > 0.001) ? y * 0.0001 : 0.000001);

                    // strike derivative
                    strikep = strike * Math.Exp(dy);
                    strikem = strike / Math.Exp(dy);
                    w = impVol.Evaluate(t, strike);
                    wp = impVol.Evaluate(t, strikep);
                    wm = impVol.Evaluate(t, strikem);
                    dwdy = (wp - wm) / (2.0 * dy);
                    d2wdy2 = (wp - 2.0 * w + wm) / (dy * dy);

                    // time derivative
                    if (t == 0.0)
                    {
                        dt = 0.0001;
                        strikept = strike * Math.Exp(integrate.AdaptLobatto(0.0, t + dt));
                        wpt = impVol.Evaluate(t + dt, strikept);
                        // if (wpt < w)
                        //    Console.WriteLine("Decreasing variance at strike {0} between time {1} and time {2}", strike, t, t + dt);
                        dwdt = (wpt - w) / dt;
                    }
                    else
                    {
                        dt = Math.Min(0.0001, t / 2.0);
                        strikept = strike * Math.Exp(integrate.AdaptLobatto(t, t + dt));
                        strikemt = strike * Math.Exp(-integrate.AdaptLobatto(t - dt, t));
                        wpt = impVol.Evaluate(t + dt, strikept);
                        wmt = impVol.Evaluate(t + dt, strikemt);

                        //if (wpt < w) 
                        //    Console.WriteLine("Decreasing variance at strike {0} between time {1} and time {2}", strike, t, t + dt);
                        //if (w < wmt) 
                        //    Console.WriteLine("Decreasing variance at strike {0} between time {1} and time {2}", strike, t-dt, t);
                        dwdt = (wpt - wmt) / (2.0 * dt);
                    }
                    if (dwdy == 0.0 && d2wdy2 == 0.0)
                        locVolMatrix[i, j] = Math.Sqrt(dwdt);
                    else
                    {
                        den1 = 1.0 - y / w * dwdy;
                        den2 = 0.25 * (-0.25 - 1.0 / w + y * y / w / w) * dwdy * dwdy;
                        den3 = 0.5 * d2wdy2;
                        locVolMatrix[i, j] = dwdt / (den1 + den2 + den3);
                        //if (locVolMatrix[i,j] < 0.0)
                        //    Console.WriteLine("Negative local vol^2 at strike {0} and time {1}; " +
                        //        "Black vol surface is not smooth enought.", strike, t);
                    }
                }
            }

            // Create dupire outputs.
            PFunction2D.PFunction2D localVol = new PFunction2D.PFunction2D(locVolMat, locVolStr, locVolMatrix);
            localVol.Parse(null);
            string[] names = new string[] { "S0" };
            Vector param = new Vector(1);
            param[0] = Hdataset.S0;
            EstimationResult result = new EstimationResult(names, param);
            //result.Objects = new object[3];
            result.Objects = new object[4];
            result.Objects[0] = this.r;
            result.Objects[1] = this.q;
            result.Objects[2] = localVol;
            result.Objects[3] = impVol;
            //Console.WriteLine("r = " + HCalData.Rate.ToString());
            //Console.WriteLine("q = " + HCalData.DividendYield.ToString());
            return result;
        }
    }
}
