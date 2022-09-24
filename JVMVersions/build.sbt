scalaVersion := "3.1.3"

fork := true

libraryDependencies += "org.scalafx" %% "scalafx" % "18.0.1-R28"

ThisBuild / assemblyMergeStrategy := { case _ => MergeStrategy.first } // {
//   case PathList("javax", "servlet", xs @ _*)         => MergeStrategy.first
//   case PathList(ps @ _*) if ps.last endsWith ".html" => MergeStrategy.first
//   case "application.conf"                            => MergeStrategy.concat
//   case "unwanted.txt"                                => MergeStrategy.discard
//   case x =>
//     val oldStrategy = (ThisBuild / assemblyMergeStrategy).value
//     oldStrategy(x)
// }